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

import java.util.concurrent.BlockingQueue;
import javax.swing.JOptionPane;

/**
 *
 * @author rob
 */
public class MSThread extends Thread{
	BlockingQueue<String> queue;
	Digester digester;
	double maxCentroids;
	static int maxQueueSize;
	static int proteinCount=0;
	public boolean writing = false;
	public boolean full = false;
	public MSThread(BlockingQueue<String> q, Digester _digester){
		queue = q;
		digester = _digester;
		maxCentroids = (Runtime.getRuntime().maxMemory() / 100.0) / MassSpec.numCpus; //in bytes
		maxQueueSize = queue.size();
	}
	
	@Override
	public void run() {
		// create one mass spec object per thread
		MassSpec massSpec = new MassSpec();
		while (!queue.isEmpty()) {
			//
			if (massSpec.getTotalCentroids() > maxCentroids){ // match RAM usage to allocated VM
				// to save RAM, write all current centroids to disk and reset
				// output scans data structure 
				full = true;
			}
			if (writing){
				massSpec.writeTempFile();
				full = false;
				writing = false;
			}
			if (!full){ // have to wait for writeout if full
				try {
					String data = queue.take();
					String[] parts = data.split("_");
					massSpec.processPeptides(digester.processProtein(parts[0], Integer.parseInt(parts[1])), Double.parseDouble(parts[2]), Integer.parseInt(parts[3]));
					proteinCount++;
					simulatorGUI.progressMonitor.setNote("Simulating protein "+proteinCount + " of " + maxQueueSize);
					simulatorGUI.progressMonitor.setProgress((int) (((double) proteinCount/(double) maxQueueSize) * 75.0)+4);
				} catch (InterruptedException e) {
					JOptionPane.showMessageDialog(null, "Error: Mass spec thread interrupted.", "Error", JOptionPane.ERROR_MESSAGE);
					return;
				}
			} else{
				try{Thread.sleep(10000);}
				catch(Exception e){}
			}
		}
		massSpec.writeTempFile();
		massSpec = null;
		System.gc();
	}
}
