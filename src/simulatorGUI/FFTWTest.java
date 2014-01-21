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
