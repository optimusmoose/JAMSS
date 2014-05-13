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
import java.io.FilenameFilter;
import java.io.IOException;
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
	boolean inRange;
	
	public static void main(String [] args){
	/*	RandomFactory.setSeed();
		// test write-in write-out
		ArrayList<Double> _isotopeIntensities = new ArrayList<Double>();
		_isotopeIntensities.add(123.45);
		_isotopeIntensities.add(678.9);
		ArrayList<Double> _isotopeMasses = new ArrayList<Double>();
		_isotopeMasses.add(123.45);
		_isotopeMasses.add(987.6);
		double _predictedIntensity = 6.78;
		double _predictedRt = 234.234623456;
		int _traceLength = 5;
		int _rtFloor = 6;
		int _rtCeil = 9;
		double _charge = 234.2346;
		int _peptideID = 7;
		int _proteinID = 1;
		String _sequence = "ASDFPOIQWER";
		int[] _charges = {3,4,5};
		IsotopicEnvelope ie = new IsotopicEnvelope(_isotopeIntensities, _isotopeMasses, _predictedIntensity, _predictedRt, _traceLength, _rtFloor, _rtCeil,
		_charge, _peptideID, _proteinID, _sequence, _charges);
		
		String path = MassSpec.class.getProtectionDomain().getCodeSource().getLocation().getPath();
		String pathToClass = "";
		try {
			pathToClass = URLDecoder.decode(path, "UTF-8").replace("JAMSS.jar", "");
		} catch (UnsupportedEncodingException ex) {
			JOptionPane.showMessageDialog(null, "Error: encoding error when finding path to JAR file.", "Error", JOptionPane.ERROR_MESSAGE);
		}
		ie.serialize(pathToClass,0);
		LinkedList<IsotopicEnvelope> iesFromFile = IsotopicEnvelope.getIEsFromFile(pathToClass, _rtFloor, localRandomFactory);
		for(IsotopicEnvelope ieFromFile : iesFromFile){
			System.out.println(ieFromFile.isotopeIntensities.get(0) + " 123.45");
			System.out.println(ieFromFile.isotopeIntensities.get(1) + " 678.9");
			System.out.println(ieFromFile.isotopeMasses.get(0) + " 123.45");			
			System.out.println(ieFromFile.isotopeMasses.get(1) + " 987.6");			
			System.out.println(ieFromFile.predictedIntensity + " 6.78");			
			System.out.println(ieFromFile.predictedRt + " 234.234623456");			
			System.out.println(ieFromFile.traceLength + " 5");						
			System.out.println(ieFromFile.rtFloor + " 6");									
			System.out.println(ieFromFile.rtCeil + " 9");									
			System.out.println(ieFromFile.charge + " 234.2346");									
			System.out.println(ieFromFile.peptideID + " 7");									
			System.out.println(ieFromFile.proteinID + " 1");												
			System.out.println(ieFromFile.sequence + " ASDFPOIQWER");
			System.out.println(ieFromFile.charges[0] + " 3");
			System.out.println(ieFromFile.charges[1] + " 4");
			System.out.println(ieFromFile.charges[2] + " 5");
		} */
	}
	public static LinkedList<IsotopicEnvelope> getIEsFromFile(String pathToClass, int rtIdx, XORShiftRandom  rand){
		LinkedList<IsotopicEnvelope> ieList = new LinkedList<IsotopicEnvelope>();		

		
		// get all .ser files with this rtIdx
		File directory = new File(pathToClass + "JAMSSfiles");
		String[] myFiles = directory.list(new FilenameFilter() {
			@Override
			public boolean accept(File directory, String fileName) {
				return fileName.endsWith(".ser");
			}
		});
			
		for (String fileName : myFiles){
			if (fileName.split("_")[0].equals(Integer.toString(rtIdx))){
				// open file with this RT and read all IEs from it
				fileName = pathToClass + "JAMSSfiles" + File.separator + fileName;
				Path path = Paths.get(fileName);
				File f = path.toFile();
				if(f.exists()){
					try {
						byte[] ba = Files.readAllBytes(path);
						ByteBuffer ieCollection = ByteBuffer.wrap(ba);
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
							String _sequence = new String(seqCharArray);

							// int charge length
							int[] _charges = new int[ieCollection.getInt()];

							// int charges
							for(int i=0; i<_charges.length; i++){
								_charges[i] = ieCollection.getInt();
							}

							IsotopicEnvelope newIE = new IsotopicEnvelope(_isotopeIntensities, 
																			_isotopeMasses, 
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
								newIE.isotopeTraces[i] = new IsotopeTrace(newIE.predictedIntensity, newIE.isotopeIntensities.get(i), newIE.traceLength, newIE.predictedRt, newIE.isotopeMasses.get(i),rand.nextDouble());
//								if (newIE.isotopeTraces[i].isotopeTraceMass > MassSpec.maxMZ){MassSpec.maxMZ = newIE.isotopeTraces[i].isotopeTraceMass;}
							}	
							ieList.add(newIE);
						}
					} catch (IOException ex) {
						JOptionPane.showMessageDialog(null, "Error de-serializing isotopic envelopes.", "Error", JOptionPane.ERROR_MESSAGE);
					}
				}
			}
		}
		if(ieList.size() > 0){return ieList;}
		return null;
	}
	
	public IsotopicEnvelope(ArrayList<Double> _isotopeIntensities, 
							ArrayList<Double> _isotopeMasses, 
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
		inRange = true;
	}
	public boolean isInRange(int rt){

		if(MassSpec.oneD || (rt >= rtFloor && rt <= rtCeil)){
			return true;
		}
		return false;
	}
	public Centroid[] getIEAtRT(double rt, RandomFactory localRandomFactory, LinkedList<MS2> ms2s){
		Centroid[] result = new Centroid[isotopeTraces.length];
		for (int i=0; i<isotopeTraces.length; i++){
			if(isotopeTraces[i].isotopeTraceMass > MassSpec.minMZ && isotopeTraces[i].isotopeTraceMass < MassSpec.maxMZ){
				if (MassSpec.oneD && localRandomFactory.localRand.nextBoolean()) {result[i] = null;} // random drop for 1-d
				Centroid tempCent = isotopeTraces[i].getCentroidAtRT(rt, localRandomFactory);
				if (tempCent.abundance > MassSpec.minWhiteNoiseIntensity){
					tempCent.charge = (int) charge;
					tempCent.pepID = peptideID;
					tempCent.proteinID = proteinID;
					tempCent.isotopeTraceID = i;
					result[i] = tempCent;

					ms2s.add(new MS2(sequence, (int) charge, isotopeTraces[i].isotopeTraceMass, tempCent.abundance));
				} 
			}
		}
		return result;
	}
	public void serialize(String pathToClass, int id){
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
		for(int i=0; i<sequence.length(); i++){
			bb.putChar(sequence.charAt(i));
		}

		// int charge length
		bb.putInt(charges.length);
		
		// int charges
		for(int i=0; i<charges.length; i++){
			bb.putInt(charges[i]);
		}

		try{
			FileOutputStream f = new FileOutputStream(pathToClass + "JAMSSfiles" + File.separator + rtFloor + "_" + id + ".ser", true);
			bb.flip();
			byte[] data = new byte[bb.limit()];
			bb.get(data);
			f.write(data);
			f.close();
		}catch(IOException e){e.printStackTrace();}
	}
}
