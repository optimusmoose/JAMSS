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
