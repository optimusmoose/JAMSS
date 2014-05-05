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
public class IEGeneratorThread extends Thread{
	ArrayList<String> queue;
	Digester digester;
	static int totalThreads = 0;
	private int threadId;
	static int maxQueueSize;
	static int proteinCount=0;
	public boolean finished = false;
	public MassSpec massSpec;
	static public int maxEnvelopes;
	public IEGeneratorThread(ArrayList<String> q, Digester _digester, int threadIdx){
		threadId = totalThreads++;
		queue = q;
		digester = _digester;
//		if (maxQueueSize < queue.size()){
//			maxQueueSize = queue.size();
//		}
		massSpec = new MassSpec(threadIdx);
		maxEnvelopes = (int) ((Runtime.getRuntime().maxMemory() - 1000000000) / 1000) / MassSpec.numCpus;
	}
	
	
	@Override
	public void run() {
		// create one mass spec object per thread
		for(String data:queue) {
			if(Runtime.getRuntime().freeMemory() < Runtime.getRuntime().maxMemory()*0.10){
				massSpec.writeEnvelopes(threadId);
			}
			if(data==null){return;}
			String[] parts = data.split("_");
			massSpec.processPeptides(digester.processProtein(parts[0], Integer.parseInt(parts[1])), Double.parseDouble(parts[2]), Integer.parseInt(parts[3]));
			proteinCount++;
			simulatorGUI.progressMonitor.setNote("Simulating protein "+proteinCount + " of " + maxQueueSize);
			simulatorGUI.progressMonitor.setProgress((int) (((double) proteinCount/(double) maxQueueSize) * 21.0)+4);
		}
		massSpec.writeEnvelopes(threadId);
		finished = true;
		return;
	}
}
