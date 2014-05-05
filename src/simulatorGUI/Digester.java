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
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.nio.charset.Charset;
import java.util.ArrayList;
import javax.swing.JOptionPane;

/**
 *
 * @author rob
 */
public class Digester {

	private String name;
	private String splitLetters;
	private String excludeLetters;
	private Boolean cTerm;
	
	public Digester(String _name, String _splitLetters, String _excludeLetters, Boolean _cTerm) {
		name = _name;
		splitLetters = _splitLetters;
		excludeLetters = _excludeLetters;
		cTerm = _cTerm;
	}

	public String getName() {
		return name;
	}

	public String getSplitLetters() {
		return splitLetters;
	}

	public String getExcludeLetters() {
		return excludeLetters;
	}

	public Boolean getcTerm() {
		return cTerm;
	}

	public void processFile(File file, int missedCleavages, String intensityModelLocation, String rtModelLocation) {
		simulatorGUI.progressMonitor.setNote("Reading .fasta file");
		ArrayList<String> queue = new ArrayList<String>();
		
		// set up the random seed
		RandomFactory.setSeed();
				
		Charset encoding = Charset.defaultCharset();
		Reader reader;
		int fastaCTR = 0;
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
							System.out.println("Processed line " + fastaCTR + " of fasta");
							fastaCTR++;
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
		
		// one thread per core
		IEGeneratorThread[] threads = new IEGeneratorThread[MassSpec.numCpus];
		simulatorGUI.progressMonitor.setNote("Digesting and running sample through mass spectrometer.");
		simulatorGUI.progressMonitor.setProgress(2);
		
		int to = 0;
		int from = 0;
		boolean finished = false;
		int chunk = (queue.size() / MassSpec.numCpus)/2;
		IEGeneratorThread.maxQueueSize = queue.size();
		while(!finished){
			finished = true;
			for(int i = 0; i < MassSpec.numCpus; i++){
				if(threads[i] == null){
					if(to < queue.size()-1){
						finished = false;
						to = Math.min(queue.size()-1,from + chunk);
						threads[i] = new IEGeneratorThread(new ArrayList<String>(queue.subList(from,to+1)), this, i); //must do to+1 because to is exclusive
						threads[i].start();
						from = to;
					}
				} else{
					if (threads[i].finished){
						for(IsotopicEnvelope ie : threads[i].massSpec.isotopicEnvelopesInstance){
							MassSpec.isotopicEnvelopes.add(ie);
						}
						threads[i] = null;
						System.gc();
					} else{
						finished = false;
					}
				}
			}
			try{Thread.sleep(5000);} catch(Exception e){
				JOptionPane.showMessageDialog(null, "Error finishing mass spec.", "Error", JOptionPane.ERROR_MESSAGE);
			}
		}
		queue = null;
		System.gc();
		MassSpec.finishController();
	}

	public ArrayList<String> processProtein(String protein, int missedCleavages) {
		StringBuilder nextPep = new StringBuilder();
		ArrayList<String> peptides = new ArrayList<>(); // canonical
		ArrayList<String> augmentedPeptides = new ArrayList<>(); //includes missedCleavage possibilities
		for (int i = 0; i < protein.length(); i++) {
			char c1 = protein.charAt(i);
			// exclusion is true if next char exists and matches something on the ignore list
			Boolean exclusion = (i + 1 < protein.length() ? getExcludeLetters().indexOf(protein.charAt(i + 1)) != -1 : false);
			if (getSplitLetters().indexOf(c1) != -1 && !exclusion) {
				if (!getcTerm()) { //n terminus, put current char in next peptide
					if (nextPep.length() > 0) {
						peptides.add(nextPep.toString());
					}
					nextPep = new StringBuilder();
					nextPep.append(c1);
				} else { // c terminus, put current char in current peptide
					nextPep.append(c1);
					if (nextPep.length() > 0) {
						peptides.add(nextPep.toString());
					}
					nextPep = new StringBuilder();
				}
			} else {
				nextPep.append(protein.charAt(i));
			}
		}
		if (nextPep.length() > 0) {
			peptides.add(nextPep.toString());
		}
		// handle missed cleavages
		for (int i = 0; i < peptides.size(); i++) { // for each canonical peptide
			for (int j = 0; j <= missedCleavages && j <= peptides.size(); j++) { // for each possible missed cleavage (0 to missedCleavages where missedCleavages is the max)
				StringBuilder sequence = new StringBuilder();
				for (int k = i; k <= j + i && j + i < peptides.size(); k++) { // construct the series of j consecutive peptides as one sequence
					sequence.append(peptides.get(k));
				}
				if (sequence.length() > 0) {
					augmentedPeptides.add(sequence.toString());
				} // add the constructed sequence to the list of possible sequences
			}
		}
		return augmentedPeptides;
	}

	public void testProcessProtein() {
        // NOTE only for Trypsin
		// standard cleavage
		StringBuilder test = new StringBuilder();
		ArrayList<String> expected = new ArrayList<>();
		ArrayList<String> output;

		test.append("");
		if (processProtein(test.toString(), 0).size() != 0) {
			System.out.println("fail 0");
		}
		System.out.println("Done 0");

		test = new StringBuilder();
		test.append("A");
		expected = new ArrayList<>();
		expected.add("A");
		output = processProtein(test.toString(), 0);
		if (!(output.equals(expected))) {
			System.out.println("fail 1");
			System.out.println(output);
		}
		System.out.println("Done 1");

		test = new StringBuilder();
		test.append("R");
		expected = new ArrayList<>();
		expected.add("R");
		if (!(processProtein(test.toString(), 0).equals(expected))) {
			System.out.println("fail 2");
		}
		System.out.println("Done 2");

		test = new StringBuilder();
		test.append("AAA");
		expected = new ArrayList<>();
		expected.add("AAA");
		output = processProtein(test.toString(), 0);
		if (!(output.equals(expected))) {
			System.out.println("fail 3");
			System.out.println(output);
		}
		System.out.println("Done 3");

		test = new StringBuilder();
		test.append("RAA");
		expected = new ArrayList<>();
		expected.add("R");
		expected.add("AA");
		if (!(processProtein(test.toString(), 0).equals(expected))) {
			System.out.println("fail 4");
		}
		System.out.println("Done 4");

		test = new StringBuilder();
		test.append("ARA");
		expected = new ArrayList<>();
		expected.add("AR");
		expected.add("A");
		if (!(processProtein(test.toString(), 0).equals(expected))) {
			System.out.println("fail 5");
		}
		System.out.println("Done 5");

		test = new StringBuilder();
		test.append("AAR");
		expected = new ArrayList<>();
		expected.add("AAR");
		output = processProtein(test.toString(), 0);
		if (!(output.equals(expected))) {
			System.out.println("fail 6");
			System.out.println(output);
		}
		System.out.println("Done 6");

		test = new StringBuilder();
		test.append("RAR");
		expected = new ArrayList<>();
		expected.add("R");
		expected.add("AR");
		if (!(processProtein(test.toString(), 0).equals(expected))) {
			System.out.println("fail 7");
		}
		System.out.println("Done 7");

		test = new StringBuilder();
		test.append("RRR");
		expected = new ArrayList<>();
		expected.add("R");
		expected.add("R");
		expected.add("R");
		if (!(processProtein(test.toString(), 0).equals(expected))) {
			System.out.println("fail 8");
		}
		System.out.println("Done 8");

		test = new StringBuilder();
		test.append("AA");
		expected = new ArrayList<>();
		expected.add("AA");
		if (!(processProtein(test.toString(), 0).equals(expected))) {
			System.out.println("fail 10");
		}
		System.out.println("Done 10");

		// test missed cleavages
		test = new StringBuilder();
		test.append("A");
		expected = new ArrayList<>();
		expected.add("A");
		if (!(processProtein(test.toString(), 1).equals(expected))) {
			System.out.println("fail 11");
		}
		System.out.println("Done 11");

		test = new StringBuilder();
		test.append("R");
		expected = new ArrayList<>();
		expected.add("R");
		if (!(processProtein(test.toString(), 1).equals(expected))) {
			System.out.println("fail 12");
		}
		System.out.println("Done 12");

		test = new StringBuilder();
		test.append("AAA");
		expected = new ArrayList<>();
		expected.add("AAA");
		output = processProtein(test.toString(), 1);
		if (!(output.equals(expected))) {
			System.out.println("fail 13");
			System.out.println(output);
		}
		System.out.println("Done 13");

		test = new StringBuilder();
		test.append("RAA");
		expected = new ArrayList<>();
		expected.add("R");
		expected.add("RAA");
		expected.add("AA");
		output = processProtein(test.toString(), 1);
		if (!(output.equals(expected))) {
			System.out.println("fail 14");
			System.out.println(output);
		}
		System.out.println("Done 14");

		test = new StringBuilder();
		test.append("ARA");
		expected = new ArrayList<>();
		expected.add("AR");
		expected.add("ARA");
		expected.add("A");
		if (!(processProtein(test.toString(), 1).equals(expected))) {
			System.out.println("fail 15");
		}
		System.out.println("Done 15");

		test = new StringBuilder();
		test.append("AAR");
		expected = new ArrayList<>();
		expected.add("AAR");
		if (!(processProtein(test.toString(), 1).equals(expected))) {
			System.out.println("fail 16");
		}
		System.out.println("Done 16");

		test = new StringBuilder();
		test.append("RRA");
		expected = new ArrayList<>();
		expected.add("R");
		expected.add("RR");
		expected.add("R");
		expected.add("RA");
		expected.add("A");
		if (!(processProtein(test.toString(), 1).equals(expected))) {
			System.out.println("fail 17");
		}
		System.out.println("Done 17");

		System.out.println("DONE TEST");
	}
	
}
