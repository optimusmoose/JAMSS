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
import java.util.Collections;
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
  public static double maxIntensity = 9000000; //max protein intensity
	//double[] intensityHistogram = {0.64,0.36,0.26,0.16,0.08,0.05,0.03,0.015,0.005};
  double[] intensityHistogram = {0.40,0.10,0.09,0.08,0.07,0.06,0.05,0.02,0.01};
	double intensityHistogramChunk = (maxIntensity - MassSpec.minWhiteNoiseIntensity)/10.0;
	
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

	public void processFile(File file, String intensityModelLocation, String rtModelLocation) {
    // First read through .fasta file to count how many lines (for time estimation)
    // Also count how many #s there are. We want one per entry or none.
    simulatorGUI.progressMonitor.setNote("Counting .fasta lines, please wait");
    int numFastaLines = 0;
    boolean hashEncountered = false; // this is a flag for whether abundance is specified
    try {
        BufferedReader reader = new BufferedReader(new FileReader(file));
        String line = null;
        boolean lookingForHash = false;
        while ((line = reader.readLine()) !=null) {
            for(int i=0; i<line.length();i++){
                char charAt = line.charAt(i);
                if(charAt == '>'){
                    if(hashEncountered && lookingForHash){
                      // if here, means there is a missing # after at least one has been found
                      JOptionPane.showMessageDialog(null, "Fasta line(s) do not specify quantity after specifying quantity on at least one line. Specify qty with '#' for each entry or remove '#' symbols from fasta.", "Error", JOptionPane.ERROR_MESSAGE);
                      return;
                    }
                    numFastaLines++;
                    lookingForHash = true;
                } else if(lookingForHash == true && charAt == '#'){
                  lookingForHash = false;
                  if(!hashEncountered && numFastaLines > 1){
                    // if here, a hash was found after the first line but
                    // previous lines contained no hash
                    JOptionPane.showMessageDialog(null, "Fasta specifies quantity on at least one entry after not specifying quantity on at least one entry. Specify qty with '#' for each entry or remove '#' symbols from fasta.", "Error", JOptionPane.ERROR_MESSAGE);
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
    
		ArrayList<Protein> queue = new ArrayList<Protein>();
		
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
								abundance = maxIntensity * Double.parseDouble(abundanceBuilder.toString());
							} else { // need to sample abundance, it isn't provided
                // Histogram divides possible intensity range into 10 equal chunks
                // this code simulates a Pareto distribution, with highest intensity very
                // unlikely and lower intensities much more likely (the code below is 
                // akin to a uniform distribution, but intensityHistogram is hard
                // coded above to turn it into a Pareto-like distribution)
                double histRand = RandomFactory.rand.nextDouble();
                int histIdx = 0;
                while (histIdx < 9 && histRand <= intensityHistogram[histIdx]){histIdx++;}
                // within the probable histogram bin with in-bin variability
                abundance = MassSpec.minWhiteNoiseIntensity + ((double) (histIdx+1)) * intensityHistogramChunk + RandomFactory.rand.nextDouble() * intensityHistogramChunk;
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
							fastaIdx++;
							inHeader = true;
							if (currentInput.length() > 0) {
								queue.add(new Protein(currentInput.toString(),abundance,proteinID));
								currentInput = new StringBuilder();
								proteinID++;
							}
							currentInput.append(ch);
						} else if (!Character.isSpace(ch)) {
							currentInput.append(ch);
						}
					}
				}
				queue.add(new Protein(currentInput.toString(),abundance,proteinID));
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
		
		MassSpec.numCpus = Math.min(MassSpec.numCpus, queue.size()); // to avoid case where #fastas < #cores
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
    ArrayList<Peptide> peptides = new ArrayList<Peptide>();
    
    try{
      ArrayList<Peptide> finishedPeptides = new ArrayList<Peptide>();
      while(!finished){
        finished = true;
        for(int i = 0; i < MassSpec.numCpus; i++){
          if(digesterThreads[i] == null){
            to = from + chunk;
            if(to <= queue.size()){
              finished = false;
              digesterThreads[i] = new DigesterThread(new ArrayList<Protein>(queue.subList(from,to)), this, i); //must do to+1 because to is exclusive
              from = to;
              digesterThreads[i].start();
            }
          } else{
            if (digesterThreads[i].finished){
              
              // collect digested peptides
              for(Peptide peptide : digesterThreads[i].digestedProteins){
                finishedPeptides.add(peptide);
              }
              digesterThreads[i] = null;
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
      
      // sort peptides for consistent peptideIDs on clone, assign IDs and output to flat files
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
      
      HashMap<String, Integer> pepIdxHash = new HashMap<String, Integer>();
      int pepIdx = 0;
      // sort finishedPeptides so they are in the same order
      Collections.sort(finishedPeptides, new PeptideSequenceComparator());
      for(Peptide peptide : finishedPeptides){
        int currPepIdx = pepIdx;
        if(!pepIdxHash.containsKey(peptide.sequence)){
          pepIdxHash.put(peptide.sequence, pepIdx);
          peptideSequenceWriter.append(currPepIdx + "," + peptide.sequence + System.getProperty("line.separator"));
          peptide.peptideID = pepIdx;
          peptide.proteinID = -1; //this variable now has no meaning, as the same sequence peptide from different proteins will now be combined
          peptides.add(peptide);
          pepIdx++;
        } else{ // this sequence is already in the hash
          currPepIdx = pepIdxHash.get(peptide.sequence);
          double newAbundance = peptides.get(currPepIdx).abundance + peptide.abundance;
          peptides.set(currPepIdx, new Peptide(newAbundance, currPepIdx, peptide.sequence, peptide.proteinID));
        }
        proteinPeptideWriter.append(peptide.proteinID + "," + currPepIdx + System.getProperty("line.separator"));
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

    //calculate Isotopic patterns and finish simulation
    // create Isotopic Envelope objects for each peptide
    IEGeneratorThread[] threads = new IEGeneratorThread[MassSpec.numCpus];
    int threadIdx = -1;
    chunk = (int) (Math.ceil((double) peptides.size() / (double) MassSpec.numCpus));
    ArrayList<Peptide> currQueue = new ArrayList<Peptide>();
    for(int i = 0; i < peptides.size(); i++){ // once per peptide
      if(i % chunk == 0){ // have reached queue capacity, start new one
        if (currQueue.size() > 0){ // don't run on first iteration
          threads[threadIdx] = new IEGeneratorThread(currQueue, threadIdx);
        }
        currQueue = new ArrayList<Peptide>(); //make new queue
        threadIdx++;
      } else { //queue not full yet
        currQueue.add(peptides.get(i));
      }
    }
    if (currQueue.size() > 0){ //means last queue has not been assigned
      threads[threadIdx] = new IEGeneratorThread(currQueue, threadIdx);
    }

    // end digestion, begin simulation
    LocalProgressMonitor.beginIsotopePatternSimulation(peptides.size());
    
    //start threads
    for(int i=0; i<MassSpec.numCpus; i++){threads[i].start();}
    
    for(int i = 0; i < MassSpec.numCpus; i++){
      while(!threads[i].finished){
          try{Thread.sleep(500);} catch(Exception e){
        }
      }
    }
    threads = null;
		System.gc();
    
    // finish the writing out of files
    MassSpec.finishController();
	}
}
