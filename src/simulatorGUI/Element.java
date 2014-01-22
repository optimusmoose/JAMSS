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

/**
 *
 * @author rob
 */
public class Element {
	public int count;
	public int defaultCount;
	private int lowMassNumber;
	private int highMassNumber;
	public int[] massNumberArray;
	public double[] relativeAbundanceArray;
	private double monoIsotopicMass;
	public Element(int _count, int _lowMassNumber, int _highMassNumber, double _monoIsotopicMass){
		count = _count;
		defaultCount = _count;
		lowMassNumber = _lowMassNumber;
		highMassNumber = _highMassNumber;
		monoIsotopicMass = _monoIsotopicMass;
	}
	public void resetCount(){
		count = defaultCount;
	}
	public int getLowMassNumber(){
		return lowMassNumber;
	}
	public int getHighMassNumber(){
		return highMassNumber;
	}
	public double getMonoIsotopicmass(){
		return monoIsotopicMass;
	}
}
