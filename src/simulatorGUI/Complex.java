/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
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
