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

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.ByteBuffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.LinkedList;
import javax.swing.JOptionPane;

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
	int traceLength; 
	int rtFloor;
	int rtCeil;
	double charge;
	int peptideID; 
	int proteinID;
	String sequence;
	int[] charges;
	IsotopeTrace[] isotopeTraces;
			
	public static LinkedList<IsotopicEnvelope> getIEsFromFile(String pathToClass, int rtIdx){
		// open file with this RT and read all IEs from it
		try {
			LinkedList<IsotopicEnvelope> ieList = new LinkedList<IsotopicEnvelope>();
			Path path = Paths.get(pathToClass + "JAMSSfiles" + File.separator + rtIdx + ".ser");
			ByteBuffer ieCollection = ByteBuffer.wrap(Files.readAllBytes(path));
			while(ieCollection.hasRemaining()){
				// int number of entries in isotope intensities
				int isotopeIntensitiesSize = ieCollection.getInt();
				
				// double isotope intensities entries
				ArrayList<Double> _isotopeIntensities = new ArrayList<Double>(isotopeIntensitiesSize);
				for(int i=0; i<isotopeIntensitiesSize; i++){
					_isotopeIntensities.add(ieCollection.getDouble());
				}
				
				// int number of entries in isotope masses
				int isotopeMassesSize = ieCollection.getInt();
				
				// double isotope mass entries
				ArrayList<Double> _isotopeMasses = new ArrayList<Double>(isotopeMassesSize);
				for(int i=0; i<isotopeMassesSize; i++){
					_isotopeMasses.add(ieCollection.getDouble());
				}
				
				// double sdShareX				
				double _sdShareX = ieCollection.getDouble();
				
				// double sdEstimate;
				double _sdEstimate = ieCollection.getDouble();

				// double predictedIntensity;
				double _predictedIntensity = ieCollection.getDouble();

				// double predictedRt;
				double _predictedRt = ieCollection.getDouble();

				// double traceLength; 
				int _traceLength = ieCollection.getInt();

				// double rtFloor;
				int _rtFloor = ieCollection.getInt();

				// double rtCeil;
				int _rtCeil = ieCollection.getInt();

				// double charge;
				double _charge = ieCollection.getDouble();

				// int peptideID; 
				int _peptideID = ieCollection.getInt();

				// int proteinID;
				int _proteinID = ieCollection.getInt();

				// int sequence length
				int sequenceLength = ieCollection.getInt();
				char[] seqCharArray = new char[sequenceLength];
				
				// char sequence letters
				for(int i=0; i<sequenceLength; i++){
					seqCharArray[i] = ieCollection.getChar();
				}
				String _sequence = seqCharArray.toString();

				// int charge length
				int[] _charges = new int[ieCollection.getInt()];

				// int charges
				for(int i=0; i<_charges.length; i++){
					_charges[i] = ieCollection.getInt();
				}
				
				IsotopicEnvelope newIE = new IsotopicEnvelope(_isotopeIntensities, 
																_isotopeMasses, 
																_sdShareX, 
																_sdEstimate, 
																_predictedIntensity, 
																_predictedRt, 
																_traceLength,
																_rtFloor,
																_rtCeil, 
																_charge, 
																_peptideID, 
																_proteinID,
																 _sequence,
																_charges);
				newIE.isotopeTraces = new IsotopeTrace[newIE.isotopeIntensities.size()];
				for (int i = 0; i < newIE.isotopeIntensities.size(); i++){
					newIE.isotopeTraces[i] = new IsotopeTrace(newIE.sdShareX, newIE.sdEstimate, newIE.predictedIntensity, newIE.isotopeIntensities.get(i), newIE.traceLength, newIE.predictedRt, newIE.isotopeMasses.get(i));
				}	
				ieList.add(newIE);
				return ieList;
			}
		} catch (IOException ex) {
			JOptionPane.showMessageDialog(null, "Error de-serializing isotopic envelopes.", "Error", JOptionPane.ERROR_MESSAGE);
		}
		return null;
	}
	
	public IsotopicEnvelope(ArrayList<Double> _isotopeIntensities, 
							ArrayList<Double> _isotopeMasses, 
							double _sdShareX, 
							double _sdEstimate, 
							double _predictedIntensity, 
							double _predictedRt, 
							int _traceLength, 
							int _rtFloor, 
							int _rtCeil, 
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
	}
	public boolean isInRange(int rt){
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
	public void serialize(String pathToClass){
		// OUTPUT FORMATTING:
		int SIZEOFDOUBLE = 8;
		int SIZEOFINT = 4;
		int SIZEOFCHAR = 2;
		int bbSize = SIZEOFINT + isotopeIntensities.size() * SIZEOFDOUBLE +
					SIZEOFINT +isotopeMasses.size() * SIZEOFDOUBLE +
					SIZEOFDOUBLE * 9 + 
					SIZEOFINT * 5 +
					SIZEOFCHAR * sequence.length() + 
					SIZEOFINT;
		ByteBuffer bb = ByteBuffer.allocate(bbSize);

		// int number of entries in isotope intensities
		bb.putInt(isotopeIntensities.size());
		
		// double isotope intensities entries
		for(double isotopeIntensity : isotopeIntensities){
			bb.putDouble(isotopeIntensity);
		}
		
		// int number of entries in isotope masses
		bb.putInt(isotopeMasses.size());
		
		// double isotope mass entries
		for(double isotopeMass : isotopeMasses){
			bb.putDouble(isotopeMass);
		}

		// double sdShareX
		bb.putDouble(sdShareX);
		
		// double sdEstimate;
		bb.putDouble(sdEstimate);
		
		// double predictedIntensity;
		bb.putDouble(predictedIntensity);
		
		// double predictedRt;
		bb.putDouble(predictedRt);
		
		// int traceLength; 
		bb.putInt(traceLength);
		
		// int rtFloor;
		bb.putInt(rtFloor);
		
		// int rtCeil;
		bb.putInt(rtCeil);
		
		// double charge;
		bb.putDouble(charge);
		
		// int peptideID; 
		bb.putInt(peptideID);
		
		// int proteinID;
		bb.putInt(proteinID);
		
		// int sequence length
		bb.putInt(sequence.length());
		
		// char sequence letters
		bb.put(sequence.getBytes());
		
		// int charge length
		bb.putInt(charges.length);
		
		// int charges
		for(int chargeEntry : charges){
			bb.putInt(chargeEntry);
		}
		
		// write this IE out to the file that corresponds to it's start scan
		try{
		FileOutputStream fileOut = new FileOutputStream(pathToClass + "JAMSSfiles" + File.separator + rtCeil + ".ser",true); //append
		ObjectOutputStream out = new ObjectOutputStream(fileOut);
		out.write(bb.array());
		out.close();
		fileOut.close();
		} catch(Exception e){
			JOptionPane.showMessageDialog(null, "Error serializing isotopic envelopes.", "Error", JOptionPane.ERROR_MESSAGE);
		}
	}
}
