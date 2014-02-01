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

import java.util.LinkedList;

/**
 *
 * @author Rob.Smith
 */
public class FileWriterThread extends Thread{
	LinkedList<Centroid> masterScan;
	double rt;
	int scanIdx;
	public boolean finished = false;

	public FileWriterThread(){
		rt = -1;
	}
	
	public FileWriterThread(LinkedList<Centroid> _masterScan, double _rt, int _scanIdx){
		masterScan = _masterScan;
		rt = _rt;
		scanIdx = _scanIdx;
	}
	
	@Override
	public void run() {
		MassSpec.finish(rt, masterScan, scanIdx);
		finished = true;
		return;
	}
}
