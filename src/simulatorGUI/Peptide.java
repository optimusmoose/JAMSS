/*
 * Copyright (C) 2014 rob
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
 * @author rob
 */
public class Peptide {
  public String sequence;
  public int peptideID;
  public double abundance;
  int proteinID;
  
  public Peptide(double _peptideIntensity, int _peptideID, String _sequence, int _proteinID){ 
    proteinID = _proteinID;
		sequence = _sequence;
		peptideID = _peptideID;
    abundance = _peptideIntensity;
	}
  
  //for use when ID will be set later
  public Peptide(String _sequence, double _peptideIntensity, int _proteinID){ 
		sequence = _sequence;
		peptideID = -1;
    abundance = _peptideIntensity;
    proteinID = _proteinID;
	}
}
