/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
