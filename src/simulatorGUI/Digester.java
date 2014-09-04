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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.URISyntaxException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.HashMap;
import javax.swing.JOptionPane;

/**
 *
 * @author rob
 */
public class Digester {

	static String name;
	static String splitLetters;
	static String excludeLetters;
	static Boolean cTerm;
	
	public Digester(String _name, String _splitLetters, String _excludeLetters, Boolean _cTerm) {
		name = _name;
		splitLetters = _splitLetters;
		excludeLetters = _excludeLetters;
		cTerm = _cTerm;
	}

	public static String getName() {
		return name;
	}

	public static String getSplitLetters() {
		return splitLetters;
	}

	public static String getExcludeLetters() {
		return excludeLetters;
	}

	public static Boolean getcTerm() {
		return cTerm;
	}

	public void processFile(File file, int missedCleavages, String intensityModelLocation, String rtModelLocation) {
    // First read through .fasta file to count how many lines (for time estimation)
    // Also count how many #s there are. We want one per entry or none.
    simulatorGUI.progressMonitor.setNote("Counting .fasta lines, please wait");
    int numFastaLines = 0;
    try {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = null;
        boolean hashEncountered = false;
        boolean lookingForHash = false;
        while ((line = reader.readLine()) !=null) {
            for(int i=0; i<line.length();i++){
                char charAt = line.charAt(i);
                if(charAt == '>'){
                    if(hashEncountered && lookingForHash){
                      // if here, means there is a missing # after at least one has been found
                      JOptionPane.showMessageDialog(null, "Only some fasta entries specify qty. Specify qty with '#' for each entry or remove '#' symbols from fasta.", "Error", JOptionPane.ERROR_MESSAGE);
                      return;
                    }
                    numFastaLines++;
                    lookingForHash = true;
                } else if(lookingForHash == true && charAt == '#'){
                  lookingForHash = false;
                  if(!hashEncountered && i > 0){
                    // if here, a hash was found after the first line but
                    // previous lines contained no hash
                    JOptionPane.showMessageDialog(null, "Only some fasta entries specify qty. Specify qty with '#' for each entry or remove '#' symbols from fasta.", "Error", JOptionPane.ERROR_MESSAGE);
                    return;
                  }
                  hashEncountered = true;
                }
            }
        }
    } catch (FileNotFoundException e) {
        JOptionPane.showMessageDialog(null, "Fasta file not found.", "Error", JOptionPane.ERROR_MESSAGE);
    } catch (IOException e) {
        JOptionPane.showMessageDialog(null, "Error reading fasta file.", "Error", JOptionPane.ERROR_MESSAGE);
    }
    LocalProgressMonitor.beginDigestion(numFastaLines);
    
		ArrayList<String> queue = new ArrayList<String>();
		
		// set up the random seed
		RandomFactory.setSeed();
				
		Charset encoding = Charset.defaultCharset();
		Reader reader;
		int fastaIdx = 0;
		try {
			reader = new InputStreamReader(new FileInputStream(file), encoding);
			Reader buffer = new BufferedReader(reader);
			StringBuilder currentInput = new StringBuilder();
			StringBuilder abundanceBuilder = new StringBuilder();
			int proteinID = 0;
			double abundance = 0.0;
			int c;
			boolean inHeader = false;
			boolean inAbundance = false;
			try {
				while ((c = buffer.read()) != -1) { //until EOF, get the next character c
					char ch = (char) c; //convert the int to a char ch
					if (inHeader) {
						//continue until newline, get abundance if here
						if (ch == '\n') { //end of header
							inHeader = false;
							inAbundance = false;
							if (abundanceBuilder.length() > 0) {
								abundance = Double.parseDouble(abundanceBuilder.toString());
							}
							currentInput = new StringBuilder();
							abundanceBuilder = new StringBuilder();
						} else {
							if (ch == '#') {
								//start of abundance
								inAbundance = true;
							} else if (inAbundance) {
								abundanceBuilder.append(ch);
							} else { // keep adding to header
								currentInput.append(ch);
							}
						}

					} else {
						if (ch == '>') { // beg of next header
//							System.out.println("Processed line " + fastaIdx + " of fasta");
							fastaIdx++;
							inHeader = true;
							if (currentInput.length() > 0) {
								queue.add(currentInput.toString()+"_"+missedCleavages+"_"+abundance+"_"+proteinID);
								currentInput = new StringBuilder();
								proteinID++;
							}
							currentInput.append(ch);
						} else if (!Character.isSpace(ch)) {
							currentInput.append(ch);
						}
					}
				}
				queue.add(currentInput.toString()+"_"+missedCleavages+"_"+abundance+"_"+proteinID);
//TODO proteinID add to proteins_peptides with peptideID...        
				reader.close();
				buffer.close();
			} catch (IOException e) {
				JOptionPane.showMessageDialog(null, "Error reading fasta file.", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
		} catch (FileNotFoundException ex) {
			JOptionPane.showMessageDialog(null, "Error: Select a fasta file.", "Error", JOptionPane.ERROR_MESSAGE);
			return;
		}
		
		MassSpec.numCpus = Math.min(Runtime.getRuntime().availableProcessors(), queue.size()); // to avoid case where #fastas < #cores
		MassSpec.intensityModelLocation = intensityModelLocation;
		MassSpec.rtModelLocation = rtModelLocation;
		MassSpec.setUpRTArray();
		
    
    // Digest proteins into peptides
    DigesterThread[] digesterThreads = new DigesterThread[MassSpec.numCpus];
		simulatorGUI.progressMonitor.setNote("Digesting...");
		simulatorGUI.progressMonitor.setProgress(2);
		
		int to = 0;
		int from = 0;
		boolean finished = false;
		int chunk = (int) (Math.ceil((double) queue.size() / (double) MassSpec.numCpus));
		DigesterThread.setQueueSize(queue.size());
    HashMap<String, Integer> pepIdxHash = new HashMap<String, Integer>();
    HashMap<Integer,Double> pepAbundanceHash = new HashMap<Integer,Double>();
    try{
      try {
        MassSpec.pathToClass = MassSpec.class.getProtectionDomain().getCodeSource().getLocation().toURI().getRawPath().replace("JAMSS.jar","");
        // Java has a bug in that it gives an erroneous leading slash in windows on the above command. Workaround:
        MassSpec.pathToClass = MassSpec.pathToClass.replace("/",File.separator).substring(1); 
        FileOutputStream f = new FileOutputStream(MassSpec.pathToClass + File.separator + "test", true); //trigger exception if not on windows
  		} catch (Exception ex) {
        try {
          MassSpec.pathToClass = MassSpec.class.getProtectionDomain().getCodeSource().getLocation().toURI().getRawPath().replace("JAMSS.jar","");
        } catch (URISyntaxException ex1) {
          JOptionPane.showMessageDialog(null, "Error obtaining JAR directory.", "Error", JOptionPane.ERROR_MESSAGE);
        }
      }
      File directory = new File(MassSpec.pathToClass + "JAMSSfiles" + File.separator);
      directory.mkdir(); // create a new directory if doesn't exist
      
      //contains list of peptide sequences and IDs
			FileWriter peptideSequenceWriter = new FileWriter("JAMSSfiles" + File.separator + "peptide_sequences.csv", false);
      //contains list of peptide sequences and IDs
			FileWriter proteinPeptideWriter = new FileWriter("JAMSSfiles" + File.separator + "protein_peptides.csv", false);
      int pepIdx = 0;
      while(!finished){
        finished = true;
        for(int i = 0; i < MassSpec.numCpus; i++){
          if(digesterThreads[i] == null){
            to = from + chunk;

            if(to <= queue.size()){
              finished = false;
              digesterThreads[i] = new DigesterThread(new ArrayList<String>(queue.subList(from,to)), this, i); //must do to+1 because to is exclusive
              from = to;
              digesterThreads[i].start();
            }
          } else{
            if (digesterThreads[i].finished){
              // collect digested peptides, add peptide sequence, protein id to output file
              for(DigestedProtein digestedProtein : digesterThreads[i].digestedProteins){
                for(String peptideSequence : digestedProtein.peptideSequences){
                  int currPepIdx = pepIdx;
                  if(!pepIdxHash.containsKey(peptideSequence)){
                    pepIdxHash.put(peptideSequence, pepIdx);
                    pepAbundanceHash.put(pepIdx, digestedProtein.proteinAbundance);
                    peptideSequenceWriter.append(currPepIdx + "," + peptideSequence + System.getProperty("line.separator"));
                    pepIdx++;
                  } else{ // this sequence is already in the hash
                    currPepIdx = pepIdxHash.get(peptideSequence);
                    pepAbundanceHash.put(currPepIdx, pepAbundanceHash.get(currPepIdx) + digestedProtein.proteinAbundance);
                  }
                  proteinPeptideWriter.append(digestedProtein.proteinID + "," + currPepIdx + System.getProperty("line.separator"));
                }
                digesterThreads[i] = null;
              }
            } else{
              finished = false;
            }
          }
        }
        try{Thread.sleep(500);} catch(Exception e){
          JOptionPane.showMessageDialog(null, "Error digesting proteins.", "Error", JOptionPane.ERROR_MESSAGE);
          return;
        }
      }
      peptideSequenceWriter.flush();
			peptideSequenceWriter.close();
			proteinPeptideWriter.flush();
			proteinPeptideWriter.close();
    } catch (IOException e){
			JOptionPane.showMessageDialog(null, "Error writing output file(s).", "Error", JOptionPane.ERROR_MESSAGE);
      return;
		}
		queue = null;
		System.gc();

    //housekeeping: create array of peptide sequences in order of pepId
    String[] pepSeqArray = new String[pepIdxHash.size()];
    for(String sequence : pepIdxHash.keySet()){
      pepSeqArray[pepIdxHash.get(sequence)] = sequence;
    }
    
    //calculate Isotopic patterns and finish simulation
    // create Isotopic Envelope objects for each peptide
    IEGeneratorThread[] threads = new IEGeneratorThread[MassSpec.numCpus];
    int threadIdx = -1;
    chunk = (int) (Math.ceil((double) pepIdxHash.size() / (double) MassSpec.numCpus));
    ArrayList<Peptide> currQueue = new ArrayList<Peptide>();
    for(int i = 0; i < pepIdxHash.size(); i++){ // once per peptide
      if(i % chunk == 0){ // have reached queue capacity, start new one
        if (currQueue.size() > 0){ // don't run on first iteration
          threads[threadIdx] = new IEGeneratorThread(currQueue, threadIdx);
        }
        currQueue = new ArrayList<Peptide>(); //make new queue
        threadIdx++;
      } else { //queue not full yet
        currQueue.add(new Peptide(pepAbundanceHash.get(i), i, pepSeqArray[i]));
      }
    }
    if (currQueue.size() > 0){ //means last queue has not been started
      threads[threadIdx] = new IEGeneratorThread(currQueue, threadIdx);
    }

    // end digestion, begin simulation
    LocalProgressMonitor.beginIsotopePatternSimulation(pepIdxHash.size());
    
    //start threads
    for(int i=0; i<MassSpec.numCpus; i++){threads[i].start();}
    
    for(int i = 0; i < MassSpec.numCpus; i++){
      while(!threads[i].finished){
          try{Thread.sleep(500);} catch(Exception e){
        }
      }
    }
    threads = null;
    pepIdxHash = null;
    pepAbundanceHash = null;
    pepSeqArray = null;
		System.gc();
    
    // finish the writing out of files
    MassSpec.finishController();
	}
}
