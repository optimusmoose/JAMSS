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

import java.util.Random; 
/**
 *
 * @author rob
 */
public class Pareto {

/************************************************************************** 
 * This file is part of RandomNumberGenerator. 
 *  
 *  RandomNumberGenerator is free software: you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation, either version 3 of the License, or 
 *  (at your option) any later version. 
 *   
 *  RandomNumberGenerator is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 *   
 *  You should have received a copy of the GNU General Public License 
 *  along with RandomNumberGenerator.  If not, see //www.gnu.org/licenses/>. 
*************************************************************************/ 

	public static double pareto(Random r, double alpha, double xM) {         
		double v = r.nextDouble(); 
		while (v == 0){ 
			v = r.nextDouble(); 
		} 

		return xM / Math.pow(v, 1.0/alpha); 
	}

	public static double paretoBounded(Random r, double alpha, double L, double H) {     
		double u = r.nextDouble(); 
		while (u == 0){ 
			u = r.nextDouble(); 
		}             
		double x = -(u*Math.pow(H,alpha)-u*Math.pow(L,alpha)-Math.pow(H,alpha)) /  
						(Math.pow(H*L,alpha)); 
		return Math.pow(x, -1.0/alpha); 
	} 
	
}
