/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
	int maxCentroids;
	static int maxQueueSize;
	int index;
	public MSThread(BlockingQueue<String> q, Digester _digester){
		queue = q;
		digester = _digester;
		maxCentroids = (int) (Runtime.getRuntime().maxMemory() / 1000.0); //in b. Assume 1 kb per centroid (conservative)
		maxQueueSize = queue.size();
	}
	
	@Override
	public void run() {
		// create one mass spec object per thread
		MassSpec massSpec = new MassSpec();
		while (!queue.isEmpty()) {
			//
			if (massSpec.getTotalCentroids() > maxCentroids / MassSpec.numCpus){ // match RAM usage to allocated VM
				// to save RAM, write all current centroids to disk and reset
				// output scans data structure 
				massSpec.writeTempFile();
			}
			try {
				String data = queue.take();
				String[] parts = data.split("_");
				massSpec.processPeptides(digester.processProtein(parts[0], Integer.parseInt(parts[1])), Double.parseDouble(parts[2]), Integer.parseInt(parts[3]));
				simulatorGUI.progressMonitor.setProgress((int) (((double) queue.size()/(double) maxQueueSize) * 75.0)+4);
			} catch (InterruptedException e) {
				JOptionPane.showMessageDialog(null, "Error: Mass spec thread interrupted.", "Error", JOptionPane.ERROR_MESSAGE);
				return;
			}
		}
		massSpec.writeTempFile();
		massSpec = null;
		System.gc();
	}
}
