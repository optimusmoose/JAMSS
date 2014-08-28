/*
		JAMSS - MS simulator
    Copyright (C) 2014  Rob Smith 2robsmith@gmail.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
	
	public static void main(String [ ] args){
		Modifications testMods = new Modifications();
		testMods.powerSet = new ArrayList<>();
		testMods.modList = new ArrayList<Modification>();
		testMods.modList.add(new Modification("methionine",0.50));
		testMods.modList.add(new Modification("phosphorylation",0.50));
		testMods.modList.add(new Modification("pyroglutamate",0.50));
		testMods.buildPowerSet(testMods.modList);
		testMods.distributeIntensities();
		for(ArrayList<Modification> modList : testMods.powerSet){
			System.out.println("========");
			for (Modification mod : modList){
				System.out.println(mod.name + " " + mod.percent);
			}
		}
	}
	private void buildPowerSet(ArrayList<Modification> list){
		if(powerSet.size() == 0){ // add singles
			for(int i=0; i<list.size(); i++){
				ArrayList<Modification> temp = new ArrayList<>();
				temp.add(list.get(i));
				powerSet.add(temp);
			}
		}
		if (list.size() > 1){
			powerSet.add(list);
			for(int i=0; i<list.size(); i++)
			{
				ArrayList<Modification> temp = new ArrayList<>(list);
				temp.remove(i);
				if (temp.size() > 1){buildPowerSet(temp);}
			}
		}
	}
	public Modifications(){
		powerSet = new ArrayList<>();
		modList = new ArrayList<Modification>();
		if (spinnerOxidationMethioninePercent > 0){modList.add(new Modification("methionine",spinnerOxidationMethioninePercent));}
		if (spinnerPhosphorylationGainPercent > 0){	modList.add(new Modification("phosphorylation",spinnerPhosphorylationGainPercent));}
		if (spinnerPyroglutamateLossPercent > 0){modList.add(new Modification("pyroglutamate",spinnerPyroglutamateLossPercent));}
		buildPowerSet(modList);
		distributeIntensities();
	}
	public void distributeIntensities(){
		double[] products = new double[powerSet.size()];
		double sum_products = 0.0;
		int i = 0;
		for(ArrayList<Modification> mods : powerSet){
			double product = 1.0;
			for(Modification mod : mods){
				product*=mod.percent;
			}
			products[i] = product;
			sum_products += product;
			i++;
		}
		ArrayList<ArrayList<Modification>> newPowerSet = new ArrayList<ArrayList<Modification>>();
		i=0;
		for(ArrayList<Modification> mods : powerSet){
			ArrayList<Modification> tempPowerElement = new ArrayList<Modification>();
			double distributedIntensity = products[i] / sum_products;
			for(Modification mod : mods){
				Modification newMod = new Modification(mod.name, distributedIntensity);
				tempPowerElement.add(newMod); // have to do deep copy because you can't modify an element of an array list
			}
			newPowerSet.add(tempPowerElement);
			i++;
		}
		powerSet = newPowerSet;
	}
}
