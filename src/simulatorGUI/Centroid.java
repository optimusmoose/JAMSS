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

/**
 *
 * @author rob
 */
public class Centroid{
	public double mz;
	public double abundance;
	public static int totalCentroids=0;
	public int centroidID;
	public int isotopeTraceID;
	public int charge;
	public int pepID;
	public int proteinID;
	
	private Centroid(){}
	
	public Centroid(double _mz, double _abundance){
		centroidID = totalCentroids;
//System.out.println(totalCentroids);
		totalCentroids += 1;
		mz = _mz;
		abundance = _abundance;
	}
	
	public double getMZ(){
		return mz;
	}
	
	public Centroid clone(){
		Centroid cent = new Centroid();
		cent.abundance = this.abundance;
		cent.centroidID = this.centroidID;
		cent.charge = this.charge;
		cent.isotopeTraceID = this.isotopeTraceID;
		cent.mz = this.mz;
		cent.pepID = this.pepID;
		cent.proteinID = this.pepID;
		return cent;
	}
}
