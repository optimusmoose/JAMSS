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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import org.apache.commons.codec.binary.Base64;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.zip.Deflater;
import javax.swing.JOptionPane;
import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;
/**
 *
 * @author rob
 */
public class MassSpec {
	//
	// parameters
	//
	public static boolean oneD;
	public static double samplingRate;
	public static int runTime;
	public static double dropoutRate;
	public static double mergeThreshold;
	public static int highestNMS2;
	public static int whiteNoiseCount;
	public static double minWhiteNoiseIntensity;
	public static double maxWhiteNoiseIntensity;
	public static String truthFile;
	public static String intensityModelLocation;
	public static String rtModelLocation;
	public static int numCpus;
	public static String simOptions;
	public static boolean createTruthFile;
	public static String fastaFile;
	public static String mzmlFilename;
	public static String truthFilename;
	
	//
	// constants
	//
	private final double NEUTRON_MASS = 1.0086649156;
	private final double ABUNDANCE_THRESHOLD = 0.0001;
	
	// used for file ops
	private static String pathToClass;
	public int msIdx;
	static private int msIdxMaster = 0;
	static public boolean thisSimulationHasPoints = false;
	static int scanIdx = 0;
	
	//
	// WEKA ATTRIBUTES
	//
	// vars for RT prediction model
	private Classifier rtCls;
	private Instances rtData;
	private Attribute rtAttA;
	private Attribute rtAttR;
	private	Attribute rtAttN;
	private	Attribute rtAttD;
	private	Attribute rtAttB;
	private	Attribute rtAttC;
	private	Attribute rtAttE;
	private	Attribute rtAttQ;
	private	Attribute rtAttZ;
	private	Attribute rtAttG;
	private	Attribute rtAttH;
	private Attribute rtAttI;
	private Attribute rtAttL;
	private	Attribute rtAttK;
	private	Attribute rtAttM;
	private	Attribute rtAttF;
	private	Attribute rtAttP;
	private	Attribute rtAttS;
	private	Attribute rtAttT;
	private	Attribute rtAttW;
	private	Attribute rtAttY;
	private	Attribute rtAttV;
	private	Attribute rtAttJ;
	private	Attribute rtAttRT;

	//
	// VARS FOR COMPUTING THE CHARGE OF A PEPTIDE
	//
	int chargeY;
	int chargeC;
	int chargeK;
	int chargeH;
	int chargeR;
	int chargeD;
	int chargeE;
	int chargeU;
	int chargePolar;
	int chargeHydrophobic;
		
	// These are used for charge calculation
	// missing last values are intentional---non-existant
	private final static double[] RESIDUE_K = {2.18,8.95,10.53};
	private final static double[] RESIDUE_E = {2.19,9.67,4.25};
	private final static double[] RESIDUE_D = {1.88,9.60,3.65};
	private final static double[] RESIDUE_H = {1.82,9.17,6.00};
	private final static double[] RESIDUE_R = {2.17,9.04,12.48};
	private final static double[] RESIDUE_Q = {2.17,9.13};
	private final static double[] RESIDUE_N = {2.02,8.80};
	private final static double[] RESIDUE_C = {1.96,10.28,8.18};
	private final static double[] RESIDUE_T = {2.11,9.62};
	private final static double[] RESIDUE_S = {2.21,9.15};
	private final static double[] RESIDUE_W = {2.38,9.39};
	private final static double[] RESIDUE_Y = {2.20,9.11,10.07};
	private final static double[] RESIDUE_F = {1.83,9.13};
	private final static double[] RESIDUE_M = {2.28,9.21};
	private final static double[] RESIDUE_I = {2.36,9.68};
	private final static double[] RESIDUE_L = {2.36,9.60};
	private final static double[] RESIDUE_V = {2.32,9.62};
	private final static double[] RESIDUE_P = {1.99,10.96};
	private final static double[] RESIDUE_A = {2.34,9.69};
	private final static double[] RESIDUE_G = {2.34,9.60};
	private final static double[] RESIDUE_B = {1.95,9.20,3.65};
	private final static double[] RESIDUE_Z = {2.18,9.40,4.25};
	private final static double[] RESIDUE_X = {2.20,9.40};
	private final static double[] RESIDUE_U = {1.96,10.28,5.20};
	private static HashMap<Character,double[]> residueTable = new HashMap<>();
	
	//
	// Amino acid counts, used to build WEKA instances
	// for predictedIntensity and RT estimation
	//
	private int aaA;
	private int aaR;
	private int aaN;
	private int aaD;
	private int aaB;
	private int aaC;
	private int aaE;
	private int aaQ;
	private int aaZ;
	private int aaG;
	private int aaH;
	private int aaI;
	private int aaL;
	private int aaK;
	private int aaM;
	private int aaF;
	private int aaP;
	private int aaS;
	private int aaT;
	private int aaW;
	private int aaY;
	private int aaV;
	private int aaJ;
	
	public static double maxIntensity = 200000; 
	//double[] intensityHistogram = {0.80,0.61,0.46,0.34,0.24,0.16,0.09,0.05,0.02};
	double[] intensityHistogram = {0.64,0.36,0.26,0.16,0.08,0.05,0.03,0.015,0.005};
	double intensityHistogramChunk = (MassSpec.maxIntensity - MassSpec.minWhiteNoiseIntensity)/10.0;
	private RandomFactory localRandomFactory;
		
	FFTbase fftBase = new FFTbase();
	
	// this the hash of the centroids generated for the specific
	// thread containing this mass spec object
	public static LinkedList<IsotopicEnvelope> isotopicEnvelopes = new LinkedList<IsotopicEnvelope>();
	public LinkedList<IsotopicEnvelope> isotopicEnvelopesInstance = new LinkedList<IsotopicEnvelope>();
//	public static Object[] ieByRtIdx;
//	public Object[] ieByRtIdxInstance;
	
	// this is used to define the window for noise points
	// while allowing threaded operations
	public static double maxMZ = 2000;
	public static double minMZ = 1500;
	
	// these are used for generating MS2 distributions
	static MS2[][] masterMS2s;
	
	private static final HashMap<Character, Double> monoResidueMasses = new HashMap<>();
	// build monoResidueMasses table
	static {
		monoResidueMasses.put('A',71.037114);
		monoResidueMasses.put('R',156.101111);
		monoResidueMasses.put('N',114.042927);
		monoResidueMasses.put('D', 115.026943);
		monoResidueMasses.put('C', 103.009185);
		monoResidueMasses.put('E', 129.042593);
		monoResidueMasses.put('Q', 128.058578);
		monoResidueMasses.put('G', 57.021464);
		monoResidueMasses.put('H', 137.058912);
		monoResidueMasses.put('I', 113.084064);
		monoResidueMasses.put('L', 113.084064);
		monoResidueMasses.put('K', 128.094963);
		monoResidueMasses.put('M', 131.040485);
		monoResidueMasses.put('F', 147.068414);
		monoResidueMasses.put('P', 97.052764);
		monoResidueMasses.put('S', 87.032028);
		monoResidueMasses.put('T', 101.047679);
		monoResidueMasses.put('U', 150.95363);
		monoResidueMasses.put('W', 186.079313);
		monoResidueMasses.put('Y', 163.06332);
		monoResidueMasses.put('V', 99.068414);
		monoResidueMasses.put('*', 118.805716);
		monoResidueMasses.put('B', 172.048405);
		monoResidueMasses.put('X', 118.805716);
		monoResidueMasses.put('Z', 128.550585);
	}
	
	// b, b star (- 1.00782 - 17.02654), b not (-17.02654 - 1.00782 - 18.01056)
	static final double[] nTermIonMassDeltas = {-1.00782,-18.03436,-36.04492}; 
	// a, a star -(29.00273+17.02654), a not -(17.02654 + 29.00273+18.01056)
	static final double[] cTermIonMassDeltas = {-29.00273,-46.02927,-64.03983}; 
	private static final double NTERMMASS = 1.00782;
	private static final double CTERMMASS = 27.99491 - 10.9742;

	private static List<String> peptides = Collections.synchronizedList(new ArrayList<String>());
	private static List<Integer> proteinIDs = Collections.synchronizedList(new ArrayList<Integer>());
	
	//
	// VARS FOR COMPUTING THE ISOTOPIC DISTRIBUTION (MZS AND INTENSITIES) (NIST)
	//
	public static final Element elementO = new Element(1,16,18,15.99491461956);
	static {
		elementO.massNumberArray = new int[] {16,17,18};
		elementO.relativeAbundanceArray = new double[] {0.99759, 0.00037, 0.00204};
	}
	public static final Element elementN = new Element(0,14,15,14.0030740048);
	static {
		elementN.massNumberArray = new int[] {14,15};
		elementN.relativeAbundanceArray = new double[] {0.99635, 0.00365};
	}
	public static final Element elementC = new Element(0,12,13,12.0);
	static {
		elementC.massNumberArray = new int[] {12,13};
		elementC.relativeAbundanceArray = new double[] {0.9891, 0.0109};
	}		
	public static final Element elementH = new Element(0,1,2,1.00782503207);
	static {
		elementH.massNumberArray = new int[] {1,2};
		elementH.relativeAbundanceArray = new double[] {0.999844, 0.000156};
	}
	public static final Element elementS = new Element(0,32,36,31.972071);
	static {
		elementS.massNumberArray = new int[] {32, 33, 34, 36};
		elementS.relativeAbundanceArray = new double[] {0.9493, 0.0076, 0.0429, 0.0002};
	}
	public static final Element elementP = new Element(0,31,31,30.97376163);
	static {
		elementP.massNumberArray = new int[] {31};
		elementP.relativeAbundanceArray = new double[] {1.0};
	}
	public static final Element elementSe = new Element(0,74,82,79.9165213);
	static {
		elementSe.massNumberArray = new int[] {74, 76, 77, 78, 80, 82};
		elementSe.relativeAbundanceArray = new double[] {0.0089, 0.0937, 0.0763, 0.2377, 0.4961, 0.0873};
	}
	private static final Element[] elements = new Element[] {elementO, elementN, elementC, elementH, elementS, elementP, elementSe};
	private final String aasX = "ACDEFGIHKLMNOPQRSTUV";
	private final String aasB = "ND";
	private final String aasZ = "QE";

	//
	//create the retention times for the run
	//
	public static double[] rtArray;
	public static int[] rtArrayShifted;
	
	private static final Modifications modifications = new Modifications();
	
	public static void setUpRTArray(){
		int numScans = (int) (samplingRate * runTime);
		double specTime = 0;
		double sampleTime = 1.0 / samplingRate;
		rtArray = new double[numScans];
		for(int i = 0; i < numScans; i++){
			if(RandomFactory.rand.nextDouble() > dropoutRate){
				double tempRT = specTime - sampleTime/6.0 + RandomFactory.rand.nextDouble() * sampleTime/3.0;
				tempRT = (tempRT < 0 ? 0.0 : tempRT);
				rtArray[i] = tempRT;
			}
			specTime += 1.0/samplingRate;
		}
		Arrays.sort(rtArray);
		rtArrayShifted = new int[(int)(rtArray[rtArray.length-1] * 10.0)+2];
		//rtArrayShifted = new int[(int)(rtArray[rtArray.length-1] * 10.0)+1];
		for (int i=0; i<rtArray.length; i++){
			rtArrayShifted[(int) (rtArray[i] * 10)] = i;
		}

		// retrace array and back fill blank entries with last entry
		int lastEntry = rtArray.length-1;
		for (int i=rtArrayShifted.length-1; i>=0; i--){
			if (rtArrayShifted[i] == 0){
				rtArrayShifted[i] = lastEntry;
			} else {
				lastEntry = rtArrayShifted[i];
			}
		}
		masterMS2s = new MS2[rtArray.length][];
	}
	private double monoMZ;
	public MassSpec(int _msIdx){
		// create a unique ID for this MS (used to prevent race conditions on output files)
		msIdx = _msIdx;
		
		localRandomFactory = new RandomFactory(_msIdx);
		
		//String path = MassSpec.class.getProtectionDomain().getCodeSource().getLocation().getPath();
		File directory;
		try {
			pathToClass = MassSpec.class.getProtectionDomain().getCodeSource().getLocation().toURI().getRawPath().replace("JAMSS.jar","");
			// Java has a bug in that it gives an erroneous leading slash in windows on the above command. Workaround:
			pathToClass = pathToClass.replace("/",File.separator).substring(1); 
			FileOutputStream f = new FileOutputStream(pathToClass + File.separator + "test", true); //trigger exception if not on windows
			
		} catch (Exception ex) {
			try {
				pathToClass = MassSpec.class.getProtectionDomain().getCodeSource().getLocation().toURI().getRawPath().replace("JAMSS.jar","");
			} catch (URISyntaxException ex1) {
				JOptionPane.showMessageDialog(null, "Error obtaining JAR directory.", "Error", JOptionPane.ERROR_MESSAGE);
			}
//			JOptionPane.showMessageDialog(null, "Error: encoding error when finding path to JAR file.", "Error", JOptionPane.ERROR_MESSAGE);
		}
		// check for extant RT files and delete them
		directory = new File(pathToClass + "JAMSSfiles" + File.separator);
		directory.mkdir(); // create a new directory if doesn't exist
		// get all .ser files
		String[] myFiles = directory.list(new FilenameFilter() {
			@Override
			public boolean accept(File directory, String fileName) {
				return fileName.endsWith(".ser");
			}
		});

		if (myFiles != null){
			for (String fileName : myFiles){
				// delete file
				
				File file = new File(pathToClass + "JAMSSfiles" + File.separator + fileName);
				file.delete();
			}
		}
		
		// build residueTable
		residueTable.put('K',RESIDUE_K);
		residueTable.put('E',RESIDUE_E);
		residueTable.put('D',RESIDUE_D);
		residueTable.put('H',RESIDUE_H);
		residueTable.put('R',RESIDUE_R);
		residueTable.put('Q',RESIDUE_Q);
		residueTable.put('N',RESIDUE_N);
		residueTable.put('C',RESIDUE_C);
		residueTable.put('T',RESIDUE_T);
		residueTable.put('S',RESIDUE_S);
		residueTable.put('W',RESIDUE_W);
		residueTable.put('Y',RESIDUE_Y);
		residueTable.put('F',RESIDUE_F);
		residueTable.put('M',RESIDUE_M);
		residueTable.put('I',RESIDUE_I);
		residueTable.put('L',RESIDUE_L);
		residueTable.put('V',RESIDUE_V);
		residueTable.put('P',RESIDUE_P);
		residueTable.put('A',RESIDUE_A);
		residueTable.put('G',RESIDUE_G);
		residueTable.put('B',RESIDUE_B);
		residueTable.put('Z',RESIDUE_Z);
		residueTable.put('X',RESIDUE_X);
		residueTable.put('U',RESIDUE_U);
		
		// load the models
		try {
			//rtCls = (Classifier) weka.core.SerializationHelper.read(pathToClass + "JAMSS.jar" + File.separator + "simulatorGUI" + File.separator + rtModelLocation);
			rtCls = (Classifier) weka.core.SerializationHelper.read(this.getClass().getResourceAsStream(rtModelLocation) );
		} catch (Exception ex) {
			JOptionPane.showMessageDialog(null, "Error: Weka RT model could not be loaded.", "Error", JOptionPane.ERROR_MESSAGE);
            return;
		}

		FastVector rtAttributeList = new FastVector();
		rtAttA = new Attribute("A");
		rtAttributeList.addElement(rtAttA);
		rtAttR = new Attribute("R");
		rtAttributeList.addElement(rtAttR);
		rtAttN = new Attribute("N");
		rtAttributeList.addElement(rtAttN);
		rtAttD = new Attribute("D");
		rtAttributeList.addElement(rtAttD);
		rtAttB = new Attribute("B");
		rtAttributeList.addElement(rtAttB);
		rtAttC = new Attribute("C");
		rtAttributeList.addElement(rtAttC);
		rtAttE = new Attribute("E");
		rtAttributeList.addElement(rtAttE);
		rtAttQ = new Attribute("Q");
		rtAttributeList.addElement(rtAttQ);
		rtAttZ = new Attribute("Z");
		rtAttributeList.addElement(rtAttZ);
		rtAttG = new Attribute("G");
		rtAttributeList.addElement(rtAttG);
		rtAttH = new Attribute("H");
		rtAttributeList.addElement(rtAttH);
		rtAttI = new Attribute("I");
		rtAttributeList.addElement(rtAttI);
		rtAttL = new Attribute("L");
		rtAttributeList.addElement(rtAttL);
		rtAttK = new Attribute("K");
		rtAttributeList.addElement(rtAttK);
		rtAttM = new Attribute("M");
		rtAttributeList.addElement(rtAttM);
		rtAttF = new Attribute("F");
		rtAttributeList.addElement(rtAttF);
		rtAttP = new Attribute("P");
		rtAttributeList.addElement(rtAttP);
		rtAttS = new Attribute("S");
		rtAttributeList.addElement(rtAttS);
		rtAttT = new Attribute("T");
		rtAttributeList.addElement(rtAttT);
		rtAttW = new Attribute("W");
		rtAttributeList.addElement(rtAttW);
		rtAttY = new Attribute("Y");
		rtAttributeList.addElement(rtAttY);
		rtAttV = new Attribute("V");
		rtAttributeList.addElement(rtAttV);
		rtAttJ = new Attribute("J");
		rtAttributeList.addElement(rtAttJ);
		rtAttRT = new Attribute("rt");
		rtAttributeList.addElement(rtAttRT);
		rtData = new Instances("TestInstances",rtAttributeList,1);
	}
	
	
	private void resetCounts(){
		for (Element el : elements){
			el.resetCount();
		}
		//
		// Vars for charge calculation
		//
		chargeY = 0;
		chargeC = 0;
		chargeK = 0;
		chargeH = 0;
		chargeR = 0;
		chargeD = 0;
		chargeE = 0;
		chargeU = 0;
		chargePolar = 0;
		chargeHydrophobic = 0;
		
		//
		// WEKA AA COUNTS
		//
		aaA = 0;
		aaR = 0;
		aaN = 0;
		aaD = 0;
		aaB = 0;
		aaC = 0;
		aaE = 0;
		aaQ = 0;
		aaZ = 0;
		aaG = 0;
		aaH = 0;
		aaI = 0;
		aaL = 0;
		aaK = 0;
		aaM = 0;
		aaF = 0;
		aaP = 0;
		aaS = 0;
		aaT = 0;
		aaW = 0;
		aaY = 0;
		aaV = 0;
		aaJ = 0;
	}
	
	// The purpose of this helper is to 
	// handle the variable PTMs by spinning off
	// one computeIsotopicEnvelopes for each possible
	// combination of PTMs, based on their respective
	// percentages of occurrence given by the user
	public boolean processPeptide(String sequence, double abundance, int proteinID){
		// add pep to list of peptides
		peptides.add(sequence);
		int peptideID = peptides.size();
		proteinIDs.add(proteinID);
		
		if (modifications.powerSet.size() == 0){ // no variable mods
			if(!computeIsotopicEnvelopes(sequence, null,abundance, proteinID, peptideID)){return false;}
		} else {
			for (int i=0; i < modifications.powerSet.size(); i++){ // for each combination of mods
				// create peptide with the proper intensity
				if(!computeIsotopicEnvelopes(sequence, modifications.powerSet.get(i),abundance * ((Modification)modifications.powerSet.get(i).get(0)).percent, proteinID, peptideID)){ // get(0) because the percents are all the same.
						return false;
				}
			}
		}
		
		return true;
	}
		
	//
	// Create the isotopic distribution, yields an array
	// of mzs and an array of corresponding intensities
	// for a given sequence.
	//
	// Also figures out charge
	//
	// Has to reset all atom counts
	//	
	public boolean computeIsotopicEnvelopes(String sequence, ArrayList<Modification> mods, double abundance, int proteinID, int peptideID){	
			int chargeFloor = 0;
			int chargeCeil = 0; 
			
			resetCounts(); // this allows us to avoid re-allocating the elements for each new sequence
	
			// set modifications
			boolean modMethionine = false;
			boolean modPhosphorylation = false;
			boolean modPyroglutamate = false;
			if (mods != null){
				for (int i=0; i < mods.size(); i++){
					switch (mods.get(i).name){
						case "methionine":
							modMethionine = true;
							break;
						case "phosphorylation":
							modPhosphorylation = true;
							break;
						case "pyroglutamate":
							modPyroglutamate = true;
							break;
					}
				}
			}
			
			// 1. Modify poly amino acids to standard,
			// 2. Get atom counts from the amino acid,
			// 3. Apply post translational modifications to atom counts.
			for (char rawAA : sequence.toCharArray()){
				char aa = rawAA;
				
				// poly amino acids: "X" is for any (I exclude uncommon "U" and "O")
				if (aa == 'X'){
					// poly amino acids: "X" is for any (I exclude uncommon "U" and "O")
					aa = aasX.charAt(localRandomFactory.localRand.nextInt(aasX.length()));
				} else if (aa == 'B') { 
					// poly amino acids: "B" is "N" or "D"
					aa = aasB.charAt(localRandomFactory.localRand.nextInt(aasB.length()));
					aaB += 1;
				} else if (aa == 'Z') {
					// poly amino acids: "Z" is "Q" or "E"
					aa = aasZ.charAt(localRandomFactory.localRand.nextInt(aasZ.length()));
					aaZ += 1;
				}
				switch (aa) {
					// standard amino acids: (modifications included beneath case)
					case 'A': //=> { :C =>3, :H =>5 , :O =>1 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaA += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						elementC.count += 3;
						elementH.count += 5;
						elementO.count += 1;
						elementN.count += 1;
						break;
					case 'C': // => { :C =>3, :H =>5 , :O =>1 , :N =>1 , :S =>1 , :P =>0, :Se =>0 },
						aaC += 1;
						
						// charge calculation
						chargeC += 1;
						
						// base:
						elementC.count += 3;
						elementH.count += 5;
						elementO.count += 1;
						elementN.count += 1;
						elementS.count += 1;
						
						if (Modifications.CarbamidomethylationGain) {
							// PTM: carbamidomethylation, gain, static, C 2 H 3 N 1 O 1 S 0
							elementC.count += 2;
							elementH.count += 3;
							elementN.count += 1;
							elementO.count += 1;
						}
						break;
					case 'D': // => { :C =>4, :H =>5 , :O =>3 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaD += 1;
						
						// charge calculation
						chargeD += 1;
						
						elementC.count += 4;
						elementH.count += 5;
						elementO.count += 3;
						elementN.count += 1;
						break;
					case 'E': // => { :C =>5, :H =>7 , :O =>3 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaE += 1;
						
						// charge calculation
						chargeE += 1;
						
						elementC.count += 5;
						elementH.count += 7;
						elementO.count += 3;
						elementN.count += 1;
						break;
					case 'F': // => { :C =>9, :H =>9 , :O =>1 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaF += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						elementC.count += 9;
						elementH.count += 9;
						elementO.count += 1;
						elementN.count += 1;
						break;
					case 'G': //=> { :C =>2, :H =>3 , :O =>1 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaG += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						elementC.count += 2;
						elementH.count += 3;
						elementO.count += 1;
						elementN.count += 1;
						break;
					case 'I': // => { :C =>6, :H =>11 , :O =>1 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaI += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						elementC.count += 6;
						elementH.count += 11;
						elementO.count += 1;
						elementN.count += 1;
						break;
					case 'H': // => { :C =>6, :H =>7 , :O =>1 , :N =>3 , :S =>0 , :P =>0, :Se =>0 },
						aaH += 1;
						//charge calculation
						chargeH += 1;
						
						elementC.count += 6;
						elementH.count += 7;
						elementO.count += 1;
						elementN.count += 3;
						break;
					case 'K': // => { :C =>6, :H =>12 , :O =>1 , :N =>2 , :S =>0 , :P =>0, :Se =>0 },
						aaK += 1;
						
						// charge calculation
						chargeK += 1;
						
						elementC.count += 6;
						elementH.count += 12;
						elementO.count += 1;
						elementN.count += 2;
						break;
					case 'L': // => { :C =>6, :H =>11 , :O =>1 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaL += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						elementC.count += 6;
						elementH.count += 11;
						elementO.count += 1;
						elementN.count += 1;
						break;
					case 'M': // => { :C =>5, :H =>9 , :O =>1 , :N =>1 , :S =>1 , :P =>0, :Se =>0 },
						aaM += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						// base:
						elementC.count += 5;
						elementH.count += 9;
						elementO.count += 1;
						elementN.count += 1;
						elementS.count += 1;
						
						// PTM: M -> oxidation of methionine, gain, variable, O 1
						// variable
						if (modMethionine) {
							elementO.count += 1;
						}
						break;
					case 'N': // => { :C =>4, :H =>6 , :O =>2 , :N =>2 , :S =>0 , :P =>0, :Se =>0 },
						aaN += 1;
						
						// charge calculation
						chargePolar += 1;
						
						elementC.count += 4;
						elementH.count += 6;
						elementO.count += 2;
						elementN.count += 2;
						break;
					case 'O': // => { :C =>12, :H =>19 , :O =>2 , :N =>3 , :S =>0 , :P =>0, :Se =>0 },
						elementC.count += 12;
						elementH.count += 19;
						elementO.count += 2;
						elementN.count += 3;
						break;
					case 'P': // => { :C =>5, :H =>7 , :O =>1 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaP += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						elementC.count += 5;
						elementH.count += 7;
						elementO.count += 1;
						elementN.count += 1;
						break;
					case 'Q': // => { :C =>5, :H =>8 , :O =>2 , :N =>2 , :S =>0 , :P =>0, :Se =>0 },
						aaQ += 1;
						
						// charge calculation
						chargePolar += 1;
						
						// base:
						elementC.count += 5;
						elementH.count += 8;
						elementO.count += 2;
						elementN.count += 2;
						
						if (modPyroglutamate) {
							// Q -> pyroglutamate (or pyroglutamic acid) loss, variable, N 1 H 3
							elementN.count -= 1;
							elementH.count -= 3;
						}
						break;
					case 'R': // => { :C =>6, :H =>12 , :O =>1 , :N =>4 , :S =>0 , :P =>0, :Se =>0 },
						aaR += 1;
						
						// charge calculation
						chargeR += 1;
						
						elementC.count += 6;
						elementH.count += 12;
						elementO.count += 1;
						elementN.count += 4;
						break;
					case 'S': // => { :C =>3, :H =>5 , :O =>2 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaS += 1;
						
						// charge calculation
						chargePolar += 1;
						
						// base:
						elementC.count += 3;
						elementH.count += 5;
						elementO.count += 2;
						elementN.count += 1;
						
						if (modPhosphorylation) {
							// S,T,Y -> phosphorylation, gain, H 1 O 3 P 1
							elementH.count += 1;
							elementO.count += 3;
							elementP.count += 1;
						}
						break;
					case 'T': // => { :C =>4, :H =>7 , :O =>2 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaT += 1;
						
						// charge calculation
						chargePolar += 1;
						
						// base:
						elementC.count += 4;
						elementH.count += 7;
						elementO.count += 2;
						elementN.count += 1;
						
						if (modPhosphorylation) {
							// S,T,Y -> phosphorylation, gain, H 1 O 3 P 1
							elementH.count += 1;
							elementO.count += 3;
							elementP.count += 1;
						}
						break;
					case 'U': // => { :C =>3, :H =>5 , :O =>1 , :N =>1 , :S =>0 , :P =>0, :Se =>1 },
						// charge calculation
						chargeU += 1;
						
						elementC.count += 3;
						elementH.count += 5;
						elementO.count += 1;
						elementN.count += 1;
						elementSe.count += 1;
						break;
					case 'V': // => { :C =>5, :H =>9 , :O =>1 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaV += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						elementC.count += 5;
						elementH.count += 9;
						elementO.count += 1;
						elementN.count += 1;
						break;
					case 'W': // => { :C =>11, :H =>10 , :O =>1 , :N =>2 , :S =>0 , :P =>0, :Se =>0 },
						aaW += 1;
						
						// charge calculation
						chargeHydrophobic += 1;
						
						elementC.count += 11;
						elementH.count += 10;
						elementO.count += 1;
						elementN.count += 2;
						break;
					case 'Y': // => { :C =>9, :H =>9 , :O =>2 , :N =>1 , :S =>0 , :P =>0, :Se =>0 },
						aaY += 1;
						
						// charge calc
						chargeY += 1;
						
						// base:
						elementC.count += 9;
						elementH.count += 9;
						elementO.count += 2;
						elementN.count += 1;
						
						if (modPhosphorylation) {
							// S,T,Y -> phosphorylation, gain, H 1 O 3 P 1
							elementH.count += 1;
							elementO.count += 3;
							elementP.count += 1;
						}
						break;
					default:
						JOptionPane.showMessageDialog(null, "Error: Amino acid in fasta not recognized.", "Error", JOptionPane.ERROR_MESSAGE);
						return false;
				}
			}
			
			//
			// Figure out charge
			//
			double preCharge = 0;
			preCharge += -1.0 / (1.0 + Math.pow(10.0,residueTable.get(sequence.charAt(0))[1] - Modifications.pH));
			preCharge += -chargeD / (1.0+Math.pow(10.0, 3.65-Modifications.pH));
			preCharge += -chargeE / (1.0+Math.pow(10.0, 4.25-Modifications.pH));
			preCharge += -chargeC / (1.0+Math.pow(10.0, 8.18-Modifications.pH));
			preCharge += -chargeY / (1.0+Math.pow(10.0, 10.07-Modifications.pH));
			preCharge += 1.0 / (1.0 + Math.pow(10.0, Modifications.pH - residueTable.get(sequence.charAt(sequence.length()-1))[0]));
			preCharge += chargeH / (1.0+Math.pow(10.0, Modifications.pH-6.00));
			preCharge += chargeK / (1.0+Math.pow(10.0, Modifications.pH-10.53));
			preCharge += chargeR / (1.0+Math.pow(10.0, Modifications.pH-12.48));
			chargeFloor = (int) Math.floor(preCharge);
			chargeCeil = (int) Math.ceil(preCharge);

			int[] charges;
			if (chargeFloor == 0 || chargeFloor == chargeCeil){
				charges = new int[1];
				charges[0] = chargeCeil;
			} else {
				charges = new int[2];
				charges[0] = chargeFloor;
				charges[1] = chargeCeil;
			}
			int origHCount = elementH.count;
			for (double charge : charges){
				if (charge > 0){
					//
					// Tweak H count based on charge
					//
					elementH.count = origHCount + (int) charge;
					//
					// compute the isotopic distribution (MZs and intensities)
					//

					// calculate the isotopic slice based on the atom counts
					// get the lowNominal and highNominal values
					int lowNominal = 0;
					int highNominal = 0;
					for (Element el : elements){
						lowNominal += el.getLowMassNumber() * el.count;
						highNominal += el.getHighMassNumber() * el.count;
					}

					// get the fft of the vector of relative abundances for the isotopes of each element
					int nextPow2 = 1024;
					while (highNominal > nextPow2){
						nextPow2 *= 2;
					}
					double[] relativeAbundancesReal= new double[nextPow2];
					double[] relativeAbundancesImag= new double[nextPow2];

					double[] fftAbundancesReal = new double[nextPow2];
					double[] fftAbundancesImag = new double[nextPow2];

					double monoMass = 0;
					boolean firstGo = true;
					for (Element el : elements){
						// reset the arrays to compute the relative abundances
						double[][] relativeAbundances = el.getRelativeAbundanceFFT(nextPow2, fftBase);
						relativeAbundancesReal = relativeAbundances[0];
						relativeAbundancesImag = relativeAbundances[1];

						// convolve the frequencies of each element
						double prevReal = 0;
						for(int i=0; i<relativeAbundancesReal.length; i++){
							double[] power = complexPower(relativeAbundancesReal[i],relativeAbundancesImag[i],el.count); 
							if(firstGo){
								fftAbundancesReal[i] = power[0];
								fftAbundancesImag[i] = power[1];
							} else {
								prevReal = fftAbundancesReal[i];
								fftAbundancesReal[i] = fftAbundancesReal[i]*power[0]-fftAbundancesImag[i]*power[1];
								fftAbundancesImag[i] = prevReal*power[1]+fftAbundancesImag[i]*power[0];
							}
						}
						firstGo = false;
						monoMass += el.count * el.getMonoIsotopicmass();
					}

					fftBase.fft(fftAbundancesImag, fftAbundancesReal); // Inverse FFT
				
					double maxAbundance = 0;
					double totalAbundance = 0;
					int monoMZIndex = 0;
					for(int i=lowNominal; i<highNominal-1; i++){
						if(fftAbundancesReal[i] > 0){
						totalAbundance += fftAbundancesReal[i];
							if (fftAbundancesReal[i] > maxAbundance ) {
								monoMZIndex = i;
							}
						}
					}
					double normalizedAbundance = 0;
					double newMass = 0;
					double lastMass = (monoMass + charge)/charge - NEUTRON_MASS/charge;
					ArrayList<Double> isotopeMasses = new ArrayList<Double>(10);
					ArrayList<Double> isotopeIntensities = new ArrayList<Double>(10);
						
					// we only want entries with index between lowNominal and highNominal
					for(int i=lowNominal; i<highNominal-1; i++){
						// normalize
						normalizedAbundance = fftAbundancesReal[i] / totalAbundance;
						newMass = lastMass + NEUTRON_MASS / charge;
						lastMass = newMass;

						// keep if above threshold
						if (normalizedAbundance > ABUNDANCE_THRESHOLD) {
							isotopeMasses.add(newMass); 
							isotopeIntensities.add(normalizedAbundance);
							if (i == monoMZIndex) {monoMZ = newMass;}
						}
					}
					// calculate RT: create weka instance and run on model to get RTs
					Instance rtInstance = new Instance(rtData.numAttributes());
					rtInstance.setValue(rtAttA, aaA);
					rtInstance.setValue(rtAttR, aaR);
					rtInstance.setValue(rtAttN, aaN);
					rtInstance.setValue(rtAttD, aaD);
					rtInstance.setValue(rtAttB, aaB);
					rtInstance.setValue(rtAttC, aaC);
					rtInstance.setValue(rtAttE, aaE);
					rtInstance.setValue(rtAttQ, aaQ);
					rtInstance.setValue(rtAttZ, aaZ);
					rtInstance.setValue(rtAttG, aaG);
					rtInstance.setValue(rtAttH, aaH);
					rtInstance.setValue(rtAttI, aaI);
					rtInstance.setValue(rtAttL, aaL);
					rtInstance.setValue(rtAttK, aaK);
					rtInstance.setValue(rtAttM, aaM);
					rtInstance.setValue(rtAttF, aaF);
					rtInstance.setValue(rtAttP, aaP);
					rtInstance.setValue(rtAttS, aaS);
					rtInstance.setValue(rtAttT, aaT);
					rtInstance.setValue(rtAttW, aaW);
					rtInstance.setValue(rtAttY, aaY);
					rtInstance.setValue(rtAttV, aaV);
					rtInstance.setValue(rtAttJ, aaJ);
					rtData.add(rtInstance);
					double predictedRt=0;
					try{
						// TODO predictedRt = rtCls.classifyInstance(rtInstance);
						predictedRt = localRandomFactory.localRand.nextDouble() * MassSpec.runTime; // TODO switch back - testing if weka model updates after evalution or is deterministic
					} catch (Exception ex){
						JOptionPane.showMessageDialog(null, "WEKA error 1", "Error", JOptionPane.ERROR_MESSAGE);
						return false;
					}

					// get absolute instensity of this peptide
					double predictedIntensity = maxIntensity;
					if (abundance != 0){
						predictedIntensity = maxIntensity * (abundance/100.0);
					} else { 
						double histRand = localRandomFactory.localRand.nextDouble();
						int histIdx = 0;
						while (histIdx < 9 && intensityHistogram[histIdx] <= histRand){histIdx++;}
						// within the probable histogram bin with in-bin variability
						predictedIntensity = MassSpec.minWhiteNoiseIntensity + ((double) (histIdx+1)) * intensityHistogramChunk + localRandomFactory.localRand.nextDouble() * intensityHistogramChunk;
					}
					
					//
					// Create isotopic envelope
					//
					if (predictedRt > 0) { // check this upfront so as not to waste time
						double secondsToScans = 250.0 * (double) samplingRate * 60.0; //Convert # seconds (250 max) to # scans
						int traceLength =  (int) (localRandomFactory.localRand.nextDouble() * secondsToScans);
						
						// find first rtArray index greater than value (or rtArray.size-1 if at tail)
						int rtIndex = 0;
						int start = (int)(predictedRt * 10.0);
						if (predictedRt > rtArray[rtArray.length-1]){
							rtIndex = rtArray.length-1;
						} else {rtIndex = rtArrayShifted[start];}
							
						int rtFloor = Math.max(0,rtIndex-(int) (secondsToScans/4.0));
						int rtCeil = Math.min(rtArray.length-1,rtIndex+(int) (secondsToScans/2.0));
						if(isotopeMasses.size() > 0 && (isotopeMasses.get(0) > minMZ || isotopeMasses.get(isotopeMasses.size()-1) < maxMZ)){
							IsotopicEnvelope isotopicEnvelope = new IsotopicEnvelope(isotopeIntensities, 
																					isotopeMasses, 
																					predictedIntensity, 
																					predictedRt, 
																					traceLength, 
																					rtFloor, 
																					rtCeil, 
																					charge, 
																					peptideID, 
																					proteinID,
																					sequence,
																					charges);
							// Add this isotopicEnvelope to the list of all envelopes
							isotopicEnvelopesInstance.add(isotopicEnvelope);
						}
					}
			}
		}
		return true;
	}
	
	public boolean processPeptides(ArrayList<String> peptides, double abundance, int proteinID) {
		for (String peptide : peptides) {
			if(!processPeptide(peptide, abundance, proteinID)){return false;}
		}
		return true;
	}
	
	private double[] complexMultiply(double aReal, double aImag, double bReal, double bImag){
		double[] result = {0,0};
		result[0] = aReal*bReal-aImag*bImag;
		result[1] = aReal*bImag+aImag*bReal;
		return result;
	}

	private double[] complexPower(double real, double imag, int pow){
		// convert imaginary number from rectangular to polar
		double r = Math.sqrt(Math.pow(real,2.0) + Math.pow(imag,2.0));
		double theta = Math.atan2(imag,real);
		
		// raise to powth power
		r = Math.pow(r,(double) pow);
		theta = theta * (double) pow;
		
		// convert polar back to rectangular
		double[] result = {r*Math.cos(theta),r*Math.sin(theta)};
		return result;
	}
		
	// Here each RT scan's centroids are merged if they are within the
	// user-defined mz threshold for merging
	public static void merge(LinkedList<Centroid> masterScan){
		// sort the array list
		Collections.sort(masterScan, new CentroidMZComparator());

		double mergeEnd = mergeThreshold;
		LinkedList<Centroid> mergeGroup = new LinkedList<>();
		LinkedList<Centroid> merged = new LinkedList<Centroid>();
		for (int j=0; j < masterScan.size(); j++){ // for each point in scan
			Centroid point = masterScan.get(j);
			if (point.mz < mergeEnd){
				mergeGroup.add(point);
			} else {
				// process any accumulated merge points
				if (mergeGroup.size() > 0){
					if (mergeGroup.size() == 1){
						merged.add(mergeGroup.get(0));
					} else {
						// calculate weighted mz and summed abundance
						double mergedIntensity = 0.0;
						double numerator = 0.0;

						for(Centroid mergePoint : mergeGroup){
							numerator += mergePoint.mz * mergePoint.abundance;
							mergedIntensity += mergePoint.abundance;
						}
						Centroid mergedCent = new Centroid(numerator / mergedIntensity, mergedIntensity);
						Centroid oldCent = mergeGroup.get(mergeGroup.size()-1);
						//keep other properties from last centriod in merged group
						mergedCent.centroidID = oldCent.centroidID;
						mergedCent.charge = oldCent.charge;
						mergedCent.isotopeTraceID = oldCent.isotopeTraceID;
						mergedCent.pepID = oldCent.pepID;
						mergedCent.proteinID = oldCent.proteinID;
						merged.add(mergedCent);
					}
					mergeGroup = new LinkedList<>();
				}
				mergeGroup.add(point);
				mergeEnd = point.mz + mergeThreshold;
			}
		}
		if (mergeGroup.size() > 0){ // catch end edge case
			if (mergeGroup.size() == 1){
				merged.add(mergeGroup.get(0));
			} else {
				// calculate weighted mz and summed abundance
				double mergedIntensity = 0;
				double numerator = 0;
				for(Centroid mergePoint : mergeGroup){
					numerator += mergePoint.mz * mergePoint.abundance;
					mergedIntensity += mergePoint.abundance;
				}
				merged.add(new Centroid(numerator / mergedIntensity, mergedIntensity));
			}
		}
		masterScan = merged;
		System.gc(); // suggest gc to free memory taken by unmerged linked list
	}
	
	public static boolean outputPre(int totalScans){
		File outputFile = new File("JAMSSfiles" + File.separator + mzmlFilename);
		try{
			outputFile.delete(); // overwrite any previous truth file
			outputFile.createNewFile(); 
		} catch(IOException e){
			JOptionPane.showMessageDialog(null, "Error creating mzML file.", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		try {
			FileWriter outputWriter = new FileWriter("JAMSSfiles" + File.separator + mzmlFilename, true);
			outputWriter.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>" + System.getProperty("line.separator"));
			outputWriter.write("<mzML xmlns=\"http://psi.hupo.org/ms/mzml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\" xsi:schemaLocation=\"http://psi.hupo.org/ms/mzml http://psidev.info/files/ms/mzML/xsd/mzML1.1.0.xsd\" version=\"1.1.0\" id=\"ms1_and_ms2\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t<cvList count=\"3\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<cv id=\"MS\" fullName=\"Proteomics Standards Initiative Mass Spectrometry Ontology\" URI=\"http://psidev.cvs.sourceforge.net/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo\" version=\"3.29.0\"/>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<cv id=\"UO\" fullName=\"Unit Ontology\" URI=\"http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo\" version=\"12:10:2011\"/>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<cv id=\"IMS\" fullName=\"Imaging MS Ontology\" URI=\"http://www.maldi-msi.org/download/imzml/imagingMS.obo\" version=\"0.9.1\"/>"+ System.getProperty("line.separator"));
			outputWriter.write("\t</cvList>"+ System.getProperty("line.separator"));
			outputWriter.write("\t<fileDescription>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<fileContent>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t\t<userParam name=\"Simulated Options\" value=\"" + simOptions + "\" type=\"options\"/>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t</fileContent>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<sourceFileList count=\"1\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t\t<sourceFile id=\"sourcefile1\" name=\"" + fastaFile + "\" location=\"file://\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t\t</sourceFile>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t</sourceFileList>"+ System.getProperty("line.separator"));
			outputWriter.write("\t</fileDescription>"+ System.getProperty("line.separator"));
			outputWriter.write("\t<softwareList count=\"1\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<software id=\"JAMSS\" version=\"beta\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t</software>"+ System.getProperty("line.separator"));
			outputWriter.write("\t</softwareList>"+ System.getProperty("line.separator"));
			outputWriter.write("\t<instrumentConfigurationList count=\"2\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<instrumentConfiguration id=\"IC\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000031\" name=\"instrument model\"/>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t\t<componentList count=\"0\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t\t</componentList>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t</instrumentConfiguration>"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<instrumentConfiguration id=\"IC2\">" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000556\" name=\"LTQ Orbitrap XL\" value=\"\"/>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t<componentList count=\"3\">" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t<source order=\"1\">" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000073\" name=\"electrospray ionization\" value=\"\"/>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t</source>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t<analyzer order=\"1\">" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000083\" name=\"radial ejection linear ion trap\" value=\"\"/>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t</analyzer>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t<detector order=\"1\">" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000253\" name=\"electron multiplier\" value=\"\"/>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t\t</detector>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t\t</componentList>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t<softwareRef ref=\"JAMSS\"/>" + System.getProperty("line.separator"));
			outputWriter.write("\t\t</instrumentConfiguration>" + System.getProperty("line.separator"));
			outputWriter.write("\t</instrumentConfigurationList>"+ System.getProperty("line.separator"));
			outputWriter.write("\t<dataProcessingList count=\"1\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t<dataProcessing id=\"did_nothing\">"+ System.getProperty("line.separator"));
			outputWriter.write("\t\t</dataProcessing>"+ System.getProperty("line.separator"));		
			outputWriter.write("\t</dataProcessingList>"+ System.getProperty("line.separator"));
			outputWriter.write("\t<run id=\"simulated_run\" defaultInstrumentConfigurationRef=\"IC\">"+ System.getProperty("line.separator"));
			int spectrumCount = totalScans + totalScans * highestNMS2;
			outputWriter.write("\t\t<spectrumList count=\"" + spectrumCount + "\" defaultDataProcessingRef=\"did_nothing\">"+ System.getProperty("line.separator"));
			outputWriter.flush();
			outputWriter.close();
		} catch (IOException ex) {
			JOptionPane.showMessageDialog(null, "Error writing heading of mzML file.", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}
	
	public static boolean outputPost(){
		try {
			FileWriter mzMLWriter = new FileWriter("JAMSSfiles" + File.separator + mzmlFilename,true); // append

			// write tail
			mzMLWriter.write("\t\t</spectrumList>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t</run>"+ System.getProperty("line.separator"));
			mzMLWriter.write("</mzML>");
			mzMLWriter.flush();
			mzMLWriter.close();
		} catch (IOException ex) {
			JOptionPane.showMessageDialog(null, "Error writing end of mzML file.", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}
	// output mzml file
	public static int output(LinkedList<Centroid> masterScan, double rt, int rtIdx){
		try {
			FileWriter mzMLWriter = new FileWriter("JAMSSfiles" + File.separator + mzmlFilename,true); // append
			int precursorIdx;
			
			String encodedMZ = compressCentroids(masterScan,true);
			String encodedINT = compressCentroids(masterScan,false);
			
//			int binaryMZlength = getBinaryLength(masterScan,true);
			// write out lists in a spectrum tag
			mzMLWriter.write("\t\t\t<spectrum index=\"" + scanIdx + "\" id=\"scan="+rtIdx + "\" defaultArrayLength=\""+ masterScan.size() + "\">" + System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"1\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t<scanList count=\"1\">"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t<scan>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"" + rt + "\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t</scan>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t</scanList>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t<binaryDataArrayList count=\"2\">"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t<binaryDataArray encodedLength=\"" + encodedMZ.length() + "\">"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<binary>" + encodedMZ + "</binary>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t</binaryDataArray>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t<binaryDataArray encodedLength=\"" + encodedINT.length() + "\">"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t\t<binary>" + encodedINT + "</binary>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t\t</binaryDataArray>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t\t</binaryDataArrayList>"+ System.getProperty("line.separator"));
			mzMLWriter.write("\t\t\t</spectrum>"+ System.getProperty("line.separator"));
			mzMLWriter.flush();
			precursorIdx = scanIdx;
			scanIdx++;
			
			// put MS2 here for the N highest intensity sequences in the scan
			for (MS2 ms2 : masterMS2s[rtIdx]){
				if(ms2 != null){
					LinkedList<Centroid> ms2s = new LinkedList<>();
					char[] sequence = ms2.getSequence().toCharArray();
					// for each nterm/cterm sequence substrings:
					//		nterms: all consecutive substrings from the beginning of the sequence to the end (size - 2). E.g. [0], [0,1], [0,1,2], to [0,1,2,...,n-2]
					//		cterms: all consecutive substrings from the end of the sequence to the beginning (index 1, not 0). E.g. [1],[1,2],...,[1,2,...,n-1]
					// get mass associated with each char, generate an m/z for each nterm or cterm substring (as described above)
					double runningMass = NTERMMASS;
					for (int j=0; j<sequence.length-1; j++){ // n terms use indices from 0 to seq size-2
						runningMass += monoResidueMasses.get(sequence[j]);

						for (int charge=1; charge<=ms2.getCharge(); charge++){
							// n terms
							for (double ionMassDelta : nTermIonMassDeltas){
								ms2s.add(new Centroid(((runningMass + ionMassDelta + charge) / charge),500.0));
							}
						}
					}			
					runningMass = CTERMMASS;
					for (int j=sequence.length-1; j>0; j--){ // c terms use indices from size-1 to 1
						runningMass += monoResidueMasses.get(sequence[j]);
						for (int charge=1; charge<=ms2.getCharge(); charge++){
							// c terms
							for (double ionMassDelta : cTermIonMassDeltas){
								ms2s.add(new Centroid(((runningMass + ionMassDelta + charge) / charge),500.0)); // just use 500 for now. TODO do something better
							}
						}
					}

					// encode MS2 data
					String encodedMS2MZ = compressCentroids(ms2s,true);
					String encodedMS2INT = compressCentroids(ms2s,false);
					mzMLWriter.write("\t\t\t<spectrum index=\"" + scanIdx + "\" id=\"scan="+rtIdx + "\" defaultArrayLength=\""+ ms2s.size() + "\">" + System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000511\" name=\"ms level\" value=\"2\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000130\" name=\"positive scan\" value=\"\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000580\" name=\"MSn spectrum\" value=\"\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000127\" name=\"centroid spectrum\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t<scanList count=\"1\">"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000795\" name=\"no combination\" value=\"\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t<scan instrumentConfigurationRef=\"IC2\">"+ System.getProperty("line.separator"));
					//mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000498\" name=\"full scan\" value=\"\"/>" + System.getProperty("line.separator"));
					double scanStartTime = rt + 0.01 + RandomFactory.rand.nextDouble() * (samplingRate - 0.11);
					mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000016\" name=\"scan start time\" value=\"" + scanStartTime + "\" unitCvRef=\"UO\" unitAccession=\"UO:0000010\" unitName=\"second\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t</scan>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t</scanList>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t<precursorList count=\"1\">"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t<precursor spectrumRef=\"scan=" + precursorIdx + "\">"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<selectedIonList count=\"1\">"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t\t<selectedIon>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000042\" name=\"peak intensity\" value=\"" + ms2.getAbundance() + "\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000041\" name=\"charge state\" value=\"" + ms2.getCharge() + "\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000744\" name=\"selected ion m/z\" value=\"" + ms2.getMz() +"\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t\t</selectedIon>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t</selectedIonList>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<activation>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000133\" name=\"collision-induced dissociation\" value=\"\"/>");
					mzMLWriter.write("\t\t\t\t\t\t</activation>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t</precursor>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t</precursorList>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t<binaryDataArrayList count=\"2\">"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t<binaryDataArray encodedLength=\"" + encodedMS2MZ.length() + "\">"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000514\" name=\"m/z array\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<binary>" + encodedMS2MZ + "</binary>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t</binaryDataArray>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t<binaryDataArray encodedLength=\"" + encodedMS2INT.length() + "\">"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000523\" name=\"64-bit float\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000574\" name=\"zlib compression\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<cvParam cvRef=\"MS\" accession=\"MS:1000515\" name=\"intensity array\"/>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t\t<binary>" + encodedMS2INT + "</binary>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t\t</binaryDataArray>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t\t</binaryDataArrayList>"+ System.getProperty("line.separator"));
					mzMLWriter.write("\t\t\t</spectrum>"+ System.getProperty("line.separator"));
					scanIdx++;
				}
			}
			mzMLWriter.flush();
			mzMLWriter.close();
		} catch (IOException e) {
			JOptionPane.showMessageDialog(null, "Error writing mzML content.", "Error", JOptionPane.ERROR_MESSAGE);
			return -1;
		}
		return scanIdx;
	}
	
	static String compressCentroids(LinkedList<Centroid> centroids, boolean isMz){
		byte[] cinput = new byte[centroids.size() * 8];
		ByteBuffer buf = ByteBuffer.wrap(cinput);
		buf.order(ByteOrder.LITTLE_ENDIAN);
		if (isMz){
			for (Centroid cent : centroids){
				buf.putDouble(cent.mz);
			}
		} else {
			for (Centroid cent : centroids){
				buf.putDouble(cent.abundance);
			}			
		}
		
		byte[] input = buf.array();
		byte[] output = new byte[input.length * 2];
		Deflater compresser = new Deflater();
		compresser.setInput(input);
		compresser.finish();
		int compressedLength = compresser.deflate(output);
		compresser.end();
		byte[] compressed = new byte[compressedLength];
		for(int i = 0; i < compressedLength; i++){
			compressed[i] = output[i];
		}
		
		String decrypted = Base64.encodeBase64String(compressed);
		return decrypted;
	}
	
	static int getBinaryLength(LinkedList<Centroid> centroids, boolean isMz){
		byte[] cinput = new byte[centroids.size() * 8];
		ByteBuffer buf = ByteBuffer.wrap(cinput);
		buf.order(ByteOrder.LITTLE_ENDIAN);
		if (isMz){
			for (Centroid cent : centroids){
				buf.putDouble(cent.mz);
			}
		} else {
			for (Centroid cent : centroids){
				buf.putDouble(cent.abundance);
			}			
		}
		
		byte[] input = buf.array();
		return input.length;
	}
	
	static void createWhiteNoise(LinkedList<Centroid> masterScan){
		// create white noise
		// NOTE: must have at least one noise point per scan to avoid empty scans
		for(int j=0; j<whiteNoiseCount; j++){
			double mz = RandomFactory.rand.nextDouble() * (maxMZ);
			double abundance = minWhiteNoiseIntensity + RandomFactory.rand.nextDouble() * (maxWhiteNoiseIntensity - minWhiteNoiseIntensity);
			masterScan.add(new Centroid(mz, abundance));
		}
	}
	
	static boolean printTruthPre(){
		File peptideFile = new File("JAMSSfiles" + File.separator + "output_truth_peptides.csv");
		try{
			peptideFile.createNewFile(); // overwrite any previous truth file
		} catch(IOException e){
			JOptionPane.showMessageDialog(null, "Error creating peptide truth file.", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		
		File truthFile = new File("JAMSSfiles" + File.separator + "output_truth.csv");
		try{
			truthFile.createNewFile(); // overwrite any previous truth file
		} catch(IOException e){
			JOptionPane.showMessageDialog(null, "Error creating truth file.", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}
	
	static int charArrayHelper(char[] array, char[] buffer, int bufIndex){
		for(int i=0;i<array.length; i++){
			buffer[bufIndex] = array[i];
			bufIndex++;
		}
		return bufIndex;
	}
	static boolean printTruth(LinkedList<Centroid> masterScan, double rt){
		// print regular truth file
		try {
			// centroidID, isotopeTraceID, charge, pepID, proteinID, m/z, RT, intensity
			char[] lineEnd = System.getProperty("line.separator").toCharArray();
			char comma = ',';
			FileWriter writer = new FileWriter("JAMSSfiles" + File.separator + truthFilename ,true); // append
			BufferedWriter truthWriter = new BufferedWriter(writer, masterScan.size() * 25);
			int bufIndex = 0;
			char[] lineBuf = new char[350];
			for (Centroid cent : masterScan){
				char[] centId = Integer.toString(cent.centroidID).toCharArray();
				bufIndex = charArrayHelper(centId, lineBuf, bufIndex);
				lineBuf[bufIndex] = comma;
				bufIndex++;
				
				char[] traceId = Integer.toString(cent.isotopeTraceID).toCharArray();
				bufIndex = charArrayHelper(traceId,lineBuf,bufIndex);
				lineBuf[bufIndex] = comma;
				bufIndex++;
				
				char[] charge = Integer.toString(cent.charge).toCharArray();
				bufIndex = charArrayHelper(charge,lineBuf,bufIndex);
				lineBuf[bufIndex] = comma;
				bufIndex++;
				
				char[] pepId = Integer.toString(cent.pepID).toCharArray();
				bufIndex = charArrayHelper(pepId,lineBuf,bufIndex);
				lineBuf[bufIndex] = comma;
				bufIndex++;
				
				char[] protId = Integer.toString(cent.proteinID).toCharArray();
				bufIndex = charArrayHelper(protId,lineBuf,bufIndex);
				lineBuf[bufIndex] = comma;
				bufIndex++;
				
				char[] mz = Double.toString(cent.mz).toCharArray();
				bufIndex = charArrayHelper(mz,lineBuf,bufIndex);
				lineBuf[bufIndex] = comma;
				bufIndex++;
				
				char[] crt = Double.toString(rt).toCharArray();
				bufIndex = charArrayHelper(crt,lineBuf,bufIndex);
				lineBuf[bufIndex] = comma;
				bufIndex++;
				
				char[] abundance = Double.toString(cent.abundance).toCharArray();
				bufIndex = charArrayHelper(abundance,lineBuf,bufIndex);
				bufIndex = charArrayHelper(lineEnd,lineBuf,bufIndex);

				truthWriter.write(lineBuf, 0, bufIndex);
				//cent.centroidID.toCharArray() + comma + cent.isotopeTraceID + "," + cent.charge +"," + cent.pepID + "," + cent.proteinID + "," + cent.mz + "," + rt + "," + cent.abundance  + lineEnd;
				//truthWriter.write();
				bufIndex = 0;
			}

			truthWriter.flush();
			truthWriter.close();
 		} catch (IOException e) {
			JOptionPane.showMessageDialog(null, "Error writing to truth file.", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}

	static boolean printSequences(){
		try {
			//contains list of peptide sequences and IDs
			FileWriter peptideWriter = new FileWriter("JAMSSfiles" + File.separator + "output_truth_peptides.csv",true); // append

			for(int pepIdx = 0; pepIdx < peptides.size(); pepIdx++){
				//proteinID (pos in fasta file), peptideID, peptide sequence
				peptideWriter.append(proteinIDs.get(pepIdx) + "," + pepIdx + "," + peptides.get(pepIdx) + System.getProperty("line.separator"));
			}
			peptideWriter.flush();
			peptideWriter.close();
		} catch (IOException e){
			JOptionPane.showMessageDialog(null, "Error appending to peptide truth file.", "Error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
		return true;
	}
	
	public static void finishController(){
		
		// write the beginning of the output and truth files
		if(!outputPre(rtArray.length)){return;}
		if(createTruthFile){
			if (!printTruthPre()){return;}
		}

		FileWriterThread writeThread = new FileWriterThread();

		//create local randoms for each thread to ensure reproducibility when cloning
		RandomFactory[] localRandoms = new RandomFactory[numCpus];
		for(int j = 0; j < numCpus; j++){
			localRandoms[j] = new RandomFactory(j);
		}

		for(int i=0; i<rtArray.length; i++){
			simulatorGUI.progressMonitor.setNote("Calculating/writing RT scans: scan "+ i + " of " + rtArray.length);
			simulatorGUI.progressMonitor.setProgress((int) (((double) i/(double) rtArray.length) * 75.0)+25);
			ScanGeneratorThread.numFinished = 0; //reset each RT
			LinkedList<MS2> ms2s = new LinkedList<MS2>(); // highest intensity sequences from run
			
			// for each current RT, get rid of any out of range
			//ArrayList<Integer> iesToRemove = new ArrayList<Integer>(isotopicEnvelopes.size());
			ArrayList<Integer> iesToRemove = new ArrayList<Integer>();
			//for(int j=0; j<isotopicEnvelopes.size(); j++){
			int ieIndex = 0;
			for(IsotopicEnvelope ie : isotopicEnvelopes){
				if(!ie.isInRange(i)){
					iesToRemove.add(ieIndex);
				}
			}
			for(int j=iesToRemove.size()-1; j>=0; j--){
				isotopicEnvelopes.remove(j);	
			}
			// get ies that start at current scan from disk			
			LinkedList<IsotopicEnvelope> iesFromFile = IsotopicEnvelope.getIEsFromFile(pathToClass, i, RandomFactory.rand);
			if(iesFromFile != null){
				isotopicEnvelopes.addAll(iesFromFile);
			}
			// spin up #cpu threads to create centriods for this scan
			int chunk = Math.max(1, (isotopicEnvelopes.size() / numCpus / 2));
			ScanGeneratorThread[] threads = new ScanGeneratorThread[numCpus];
			ieIndex = 0;
			LinkedList<Centroid> masterScan = new LinkedList<Centroid>();
			while(ieIndex < isotopicEnvelopes.size()){
				for(int j = 0; j < numCpus && ieIndex < isotopicEnvelopes.size(); j++){
					if(threads[j] == null){ // thread not running
						ArrayList<Integer> ieIds = new ArrayList<Integer>(chunk);
						for(int k=ieIndex; k<isotopicEnvelopes.size() && k < ieIndex + chunk; k++){
							ieIds.add(k);
						}
						ieIndex += chunk;
						threads[j] = new ScanGeneratorThread(ieIds, rtArray[i], localRandoms[j]);
						threads[j].start();
					} else { // thread is active
						if (threads[j].finished){
							for(Centroid cent : threads[j].masterScan){
								masterScan.add(cent);
							}
							for(int k=0; k<highestNMS2 && k < threads[j].ms2s.size(); k++){
								ms2s.add(threads[j].ms2s.pop());
							}
							threads[j] = null;
						}
					}
				}
				System.gc();
				try{
					Thread.sleep(100);
				} catch(Exception e){
					JOptionPane.showMessageDialog(null, "Error writing out RT scans.", "Error", JOptionPane.ERROR_MESSAGE);
				}
			}
			// wait till they are done
			for(int j = 0; j < numCpus; j++){
				if(threads[j] != null){
					while(!threads[j].finished){
						try{
							Thread.sleep(4000);
						} catch(Exception e){
							JOptionPane.showMessageDialog(null, "Error writing out RT scans.", "Error", JOptionPane.ERROR_MESSAGE);
						}
					}
					for(Centroid cent : threads[j].masterScan){
						masterScan.add(cent);
					}
					for(int k=0; k<highestNMS2 && k < threads[j].ms2s.size(); k++){
						ms2s.add(threads[j].ms2s.pop());
					}
					threads[j] = null;
					System.gc();
				}
			}
			// process the ms2s down to only highestNMS2
			Collections.sort(ms2s, new MS2Comparator());
			Collections.reverse(ms2s); // biggest intensity first
			MS2[] ms2sArray = new MS2[highestNMS2];
			for(int k=0;k<highestNMS2 && k<ms2s.size(); k++){ms2sArray[k] = ms2s.pop();}
			masterMS2s[i] = ms2sArray;
						
			// only allow one scan to be writing at a time to maintain order
			if(writeThread.rt!=-1 && !writeThread.finished){ // will = -1 on first time through
				while(!writeThread.finished){
					try{Thread.sleep(2000);
					}catch(Exception e){}
				}
			}
			// spin a thread to write it out and move on to the next rt 
			writeThread = new FileWriterThread(masterScan, rtArray[i], i);
			writeThread.start();
			// if there isn't at least 2 GB of free RAM, wait until
			// previous writing threads finish up to make room for new
			// Scans
			while(!writeThread.finished){
				try{
					Thread.sleep(1000);
				} catch(Exception e){
					e.printStackTrace();
					System.exit(1);
				}
			}
		}
		while(!writeThread.finished){
			try{Thread.sleep(2000);} catch(Exception e){}
		}
		if(!thisSimulationHasPoints){
			JOptionPane.showMessageDialog(null, "Sorry, no peptides with above-noise intensity appear within the m/z range of this simulation.", "Finished", JOptionPane.ERROR_MESSAGE);
		}
			// write out the peptide sequences
			if(createTruthFile){printSequences();}
			// write the end of the output and truth files
			outputPost();
//		} else{
//			JOptionPane.showMessageDialog(null, "Sorry, no peptides with above-noise intensity appear within the m/z range of this simulation.", "Finished", JOptionPane.ERROR_MESSAGE);
//		}
	}
	public static void finish(double rt, LinkedList<Centroid> masterScan,int rtIdx){
		createWhiteNoise(masterScan);
		merge(masterScan);
		int status = output(masterScan, rt, rtIdx);
		if (status == -1){
			return;
		}
		if(createTruthFile){
			if(!printTruth(masterScan, rt)){
				return;
			}
		}
	}
	
	public void writeEnvelopes(int id){
		for(IsotopicEnvelope ie : isotopicEnvelopesInstance){
			ie.serialize(pathToClass,id);
		}
		isotopicEnvelopesInstance = new LinkedList<IsotopicEnvelope>();
	}
}
