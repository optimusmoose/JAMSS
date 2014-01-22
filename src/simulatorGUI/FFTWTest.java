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
public class FFTWTest {
	public static void main(String[] args){
		double[] power = complexPower(6.6,0.0,2);
		System.out.println(power[0] +" " + power[1]);
		power = complexPower(-1.65,0.952628,2);
		System.out.println(power[0] +" " + power[1]);
		power = complexPower(5.1,-0.1,5);
		System.out.println(power[0] +" " + power[1]);
		
	}
	private static double[] complexPower(double real, double imag, int pow){
		// convert imaginary number from rectangular to polar
		double r = Math.sqrt(Math.pow(real,2.0) + Math.pow(imag,2.0));
		double theta = Math.abs(Math.tanh(imag/real));
		
		// raise to powth power
		r = Math.pow(r,pow);
		theta = theta * pow;
		
		// convert polar back to rectangular
		double[] result = {r*Math.cos(theta),r*Math.sin(theta)};
		return result;
	}
}
