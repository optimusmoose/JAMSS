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

/**
 *
 * @author Rob.Smith
 */
public class MS2 {
	private String sequence;
	private int charge;
	private double mz;
	private double abundance;
	public MS2(String _sequence, int _charge, double _mz, double _abundance){
		sequence = _sequence;
		charge = _charge;
		mz = _mz;
		abundance = _abundance;
	}
	public double getAbundance(){
		return abundance;
	}
	public String getSequence(){
		return sequence;
	}
	public double getMz(){
		return mz;
	}
	public int getCharge(){
		return charge;
	}
}
