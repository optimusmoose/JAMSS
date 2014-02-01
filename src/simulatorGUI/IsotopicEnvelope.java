/*
 * Copyright (C) 2014 Rob.Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package simulatorGUI;

import java.util.ArrayList;
import java.util.LinkedList;

/**
 *
 * @author Rob.Smith
 */
public class IsotopicEnvelope {
	ArrayList<Double> isotopeIntensities; 
	ArrayList<Double> isotopeMasses;
	double sdShareX;
	double sdEstimate;
	double predictedIntensity;
	double predictedRt;
	double traceLength; 
	double rtFloor;
	double rtCeil;
	double charge;
	int peptideID; 
	int proteinID;
	String sequence;
	int[] charges;
	IsotopeTrace[] isotopeTraces;
			
	public IsotopicEnvelope(ArrayList<Double> _isotopeIntensities, 
							ArrayList<Double> _isotopeMasses, 
							double _sdShareX, 
							double _sdEstimate, 
							double _predictedIntensity, 
							double _predictedRt, 
							double _traceLength, 
							double _rtFloor, 
							double _rtCeil, 
							double _charge, 
							int _peptideID, 
							int _proteinID,
							String _sequence,
							int[] _charges){
		isotopeIntensities = _isotopeIntensities; 
		isotopeMasses = _isotopeMasses;
		sdShareX = _sdShareX;
		sdEstimate = _sdEstimate;
		predictedIntensity = _predictedIntensity;
		predictedRt = _predictedRt;
		traceLength = _traceLength; 
		rtFloor = _rtFloor;
		rtCeil = _rtCeil;
		charge = _charge;
		peptideID = _peptideID; 
		proteinID = _proteinID;
		sequence = _sequence;
		charges = _charges;
		isotopeTraces = new IsotopeTrace[isotopeIntensities.size()];
		for (int i = 0; i < isotopeIntensities.size(); i++){
			isotopeTraces[i] = new IsotopeTrace(sdShareX, sdEstimate, predictedIntensity, isotopeIntensities.get(i), traceLength, predictedRt, isotopeMasses.get(i));
		}		

	}
	public boolean isInRange(double rt){
		if(MassSpec.oneD || (rt >= rtFloor && rt <= rtCeil)){return true;}
		return false;
	}
	public Centroid[] getTraceAtRT(double rt){
		Centroid[] result = new Centroid[isotopeTraces.length];
		for (int i=0; i<isotopeTraces.length; i++){
			if (MassSpec.oneD && RandomFactory.rand.nextBoolean()) {result[i] = null;} // random drop for 1-d
			Centroid tempCent = isotopeTraces[i].getCentroidAtRT(rt);
			if (tempCent.abundance != 0){
				tempCent.charge = (int) charge;
				tempCent.pepID = peptideID;
				tempCent.proteinID = proteinID;
				tempCent.isotopeTraceID = i;
				result[i] = tempCent;

				//////////////////////////////////////////////////////
				// If this intensity is in the top N for this scan, 
				// include it in the list of sequences to be targeted
				// for MS/MS.
				if (MassSpec.highestNMS2 > 0){
					if (!MassSpec.highestNSequences.containsKey(rt)){
						// if RT not in highestNSequences, add it
						String[] sequences = new String[MassSpec.highestNMS2];
						sequences[0] = sequence;
						MassSpec.highestNSequences.put(rt, sequences);
						int[] highestCharges = new int[MassSpec.highestNMS2];
						highestCharges[0] = charges[1]; //charges[1] is the highest charge for this seq
						MassSpec.highestNCharges.put(rt, highestCharges);
						double[] intensities = new double[MassSpec.highestNMS2];
						intensities[0] = tempCent.abundance;
						MassSpec.highestNIntensities.put(rt, intensities);
						double[] mzs = new double[MassSpec.highestNMS2];
						mzs[0] = tempCent.mz;
						MassSpec.highestNMzs.put(rt, mzs);
					}
					int lowestIntensity = 0;
					//find lowest intensity in highest intensities in this scan
					for (int k=0; k<MassSpec.highestNMS2; k++){ 
						if (((double[]) MassSpec.highestNIntensities.get(rt))[k] < lowestIntensity ){
							lowestIntensity = k;
						}
					}

					// if intensity of this centroid is >= lowest intensity in highestNSequences
					if (tempCent.abundance > ((double[]) MassSpec.highestNIntensities.get(rt))[lowestIntensity]){
						((double[]) MassSpec.highestNIntensities.get(rt))[lowestIntensity] = tempCent.abundance;
						((String[]) MassSpec.highestNSequences.get(rt))[lowestIntensity] = sequence;
						((int[]) MassSpec.highestNCharges.get(rt))[lowestIntensity] = charges[1];
						((double[]) MassSpec.highestNMzs.get(rt))[lowestIntensity] = tempCent.mz;
					}
				}
			}
		}
		return result;
	}
}
