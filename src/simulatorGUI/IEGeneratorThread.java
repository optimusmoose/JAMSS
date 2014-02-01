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

/**
 *
 * @author rob
 */
public class IEGeneratorThread extends Thread{
	BlockingQueue<String> queue;
	Digester digester;
	static int maxQueueSize=0;
	static int proteinCount=0;
	public boolean finished = false;
	public MassSpec massSpec;
	public IEGeneratorThread(BlockingQueue<String> q, Digester _digester){
		queue = q;
		digester = _digester;
		if (maxQueueSize < queue.size()){
			maxQueueSize = queue.size();
		}
		massSpec = new MassSpec();
	}
	
	@Override
	public void run() {
		// create one mass spec object per thread
		while (!queue.isEmpty()) {
			String data = queue.poll();
			if(data==null){return;}
			String[] parts = data.split("_");
			massSpec.processPeptides(digester.processProtein(parts[0], Integer.parseInt(parts[1])), Double.parseDouble(parts[2]), Integer.parseInt(parts[3]));
			proteinCount++;
			simulatorGUI.progressMonitor.setNote("Simulating protein "+proteinCount + " of " + maxQueueSize + "           ");
			simulatorGUI.progressMonitor.setProgress((int) (((double) proteinCount/(double) maxQueueSize) * 40.0)+4);
		}
		finished = true;
		return;
	}
}
