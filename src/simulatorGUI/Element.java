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

import java.util.HashMap;

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
	private HashMap fftArrayHash;
	public Element(int _count, int _lowMassNumber, int _highMassNumber, double _monoIsotopicMass){
		count = _count;
		defaultCount = _count;
		lowMassNumber = _lowMassNumber;
		highMassNumber = _highMassNumber;
		monoIsotopicMass = _monoIsotopicMass;
		fftArrayHash = new HashMap();
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
	
	// Find the FFT for this element at the given length and
	// save it so we don't have to ever recompute it.
	public double[][] getRelativeAbundanceFFT(int length){
		if(!fftArrayHash.containsKey(length)){
			double[] relativeAbundancesReal= new double[length];
			double[] relativeAbundancesImag= new double[length];

			// set the relative abundance at this mass number to the element's relative abundance
			for (int i=0; i<massNumberArray.length; i++){
				relativeAbundancesReal[massNumberArray[i]] = relativeAbundanceArray[i];
			}

			// convert to frequency domain
			FFTbase fftBase = new FFTbase(relativeAbundancesReal.length);
			fftBase.fft(relativeAbundancesReal, relativeAbundancesImag);

			double[][] fftArrays = new double[2][];
			fftArrays[0] = relativeAbundancesReal;
			fftArrays[1] = relativeAbundancesImag;
			fftArrayHash.put(length, fftArrays);
		}
		double[][] fftArrays = (double[][]) fftArrayHash.get(length);
		double[][] fftArraysCopy = new double[2][length];
		for(int i=0; i<length;i++){
			fftArraysCopy[0][i] = fftArrays[0][i];
			fftArraysCopy[1][i] = fftArrays[1][i];
		}
		return fftArraysCopy;
	}
}
