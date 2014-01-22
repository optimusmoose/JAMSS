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
public class Complex {
	double real;
	double imag;
	public Complex(double _real, double _imag){
		real = _real;
		imag = _imag;
	}
	public Complex times(Complex comp){
		return new Complex(this.real * comp.real - this.imag*comp.imag, this.real*comp.imag+this.imag*comp.real);
	}
	public Complex plus(Complex comp){
		return new Complex(this.real+comp.real,this.imag+comp.imag);
	}
	public Complex minus(Complex comp){
		return new Complex(this.real-comp.real,this.imag-comp.imag);
	}
}
