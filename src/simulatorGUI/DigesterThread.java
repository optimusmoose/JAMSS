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

import java.util.ArrayList;

/**
 *
 * @author rob
 */
public class DigesterThread extends Thread{
	ArrayList<Protein> queue;
	Digester digester;
	static int maxQueueSize = 0;
	public boolean finished = false;
	static public int maxEnvelopes;
  public ArrayList<Peptide> digestedProteins;
  RandomFactory localRandomFactory;
  
	public DigesterThread(ArrayList<Protein> q, Digester _digester, int threadID){
		queue = q;
		digester = _digester;
    digestedProteins = new ArrayList<Peptide>();
		maxEnvelopes = (int) ((Runtime.getRuntime().maxMemory() - 1000000000) / 1000) / MassSpec.numCpus;
    localRandomFactory = new RandomFactory(threadID);
	}
	
	static public void setQueueSize(int size){
		if (size > maxQueueSize){
			maxQueueSize = size;
		}
	}
	
	@Override
  @SuppressWarnings("static-access")
	public void run() {
		for(Protein prot:queue) {
      digestedProteins.addAll(prot.processProtein(localRandomFactory));
      LocalProgressMonitor.updateDigestion(prot.proteinID);
		}
		finished = true;
		return;
	}
  
  /*
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
		expected = new ArrayList<>();truth
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
	} */
  
}
