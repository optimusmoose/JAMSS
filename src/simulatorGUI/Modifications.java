/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package simulatorGUI;

import java.util.ArrayList;

/**
 *
 * @author rob
 */
public class Modifications {
	public static boolean CarbamidomethylationGain;
	public static double pH;
	public static double spinnerOxidationMethioninePercent;
	public static double spinnerPhosphorylationGainPercent;
	public static double spinnerPyroglutamateLossPercent;
	public ArrayList<ArrayList<Modification>> powerSet;
	private ArrayList<Modification> modList;
	
	private void buildPowerSet(ArrayList<Modification> list, int count, double removeSum){
		if (list.size() > 0){
			powerSet.add(list);
			double percentage = 1;
			for (Modification mod : list){
				percentage *= mod.percent;
			}
			for (Modification mod : list){
				mod.percent = percentage - removeSum;
			}

			for(int i=0; i<list.size(); i++)
			{
				ArrayList<Modification> temp = new ArrayList<>(list);
				temp.remove(i);
				if (temp.size() > 0){buildPowerSet(temp, temp.size(), removeSum + percentage);}
			}
		}
	}
	public Modifications(){
		powerSet = new ArrayList<>();
		modList = new ArrayList<Modification>();
		if (spinnerOxidationMethioninePercent > 0){modList.add(new Modification("methionine",spinnerOxidationMethioninePercent));}
		if (spinnerPhosphorylationGainPercent > 0){	modList.add(new Modification("phosphorylation",spinnerPhosphorylationGainPercent));}
		if (spinnerPyroglutamateLossPercent > 0){modList.add(new Modification("pyroglutamate",spinnerPyroglutamateLossPercent));}
		buildPowerSet(modList,modList.size(), 0);
	}
}
