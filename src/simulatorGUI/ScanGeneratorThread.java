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

import java.util.ArrayList;
import java.util.LinkedList;

/**
 *
 * @author Rob.Smith
 */
public class ScanGeneratorThread extends Thread{
	ArrayList<Integer> queue;
	LinkedList<Centroid> masterScan;
	boolean finished = false;
	double rt;
	static int numFinished = 0;
	RandomFactory localRandomFactory;
	public ScanGeneratorThread(ArrayList<Integer> q, double _rt, RandomFactory _localRandomFactory){
		rt = _rt;
		queue = q;
		masterScan = new LinkedList<Centroid>();
		localRandomFactory = _localRandomFactory;
	}
	
	@Override
	public void run() {
		for(int poll : queue) {
			IsotopicEnvelope thisEnvelope = (IsotopicEnvelope) MassSpec.isotopicEnvelopes.get(poll);
			if (rt > 0){
				Centroid[] centroids = thisEnvelope.getIEAtRT(rt, localRandomFactory);
				for(int i=0; i<centroids.length; i++){
					if(centroids[i] != null){masterScan.add(centroids[i]);}
				}
				numFinished++;
			}
		}
		finished = true;
		return;
	}
	public void freeMemory(){
		masterScan = null;
	}
}
