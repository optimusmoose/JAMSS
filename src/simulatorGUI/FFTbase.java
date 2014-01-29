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
 
import java.util.HashMap;

public class FFTbase {
	int n, m;
  
  // Lookup tables.  Only need to recompute when size of FFT changes.
  double[] cos;
  double[] sin;

  static HashMap cosHash = new HashMap();
  static HashMap sinHash = new HashMap();
  
  //double[] window;
	public static void main(String [ ] args)
	{
      double[] inputReal1 = {6.4,5.2,1.3,5.3,0,0,0,0};
	  double[] inputImag1 = {0,0,0,0,0,0,0,0};
/*	  
	 double[] inputReal1 = {5,2,0,0};
	  double[] inputImag1 = {0,0,0,0};
	  double[] inputReal2 = {0,0,2,5};
	  double[] inputImag2 = {0,0,0,0};
	  
	  FFTbase me = new FFTbase(inputReal1.length);
	  me.fft(inputReal1, inputImag1);
	  me.fft(inputReal2, inputImag2);
	  for(int i=0; i<inputReal1.length; i++){
		double prevReal = inputReal1[i];
		inputReal1[i] = inputReal1[i]*inputReal2[i]-inputImag1[i]*inputImag2[i];
		inputImag1[i] = prevReal*inputImag2[i]+inputImag1[i]*inputReal2[i];
	  }
*/  
	FFTbase me = new FFTbase(inputReal1.length);
	me.fft(inputReal1, inputImag1);
	  me.fft(inputImag1,inputReal1);
	  me.printArray(inputReal1);
	  me.printArray(inputImag1);
	}
	private void printArray(double[] a){
		StringBuffer sb = new StringBuffer();
		for (int i = 0; i < a.length; i++){
			sb.append(a[i]);
			sb.append(",");
		}
		sb.deleteCharAt(sb.length()-1);
		System.out.println(sb.toString());
	}
  
	public FFTbase(int n) {
		this.n = n;
		this.m = (int)(Math.log(n) / Math.log(2));
		
		// Make sure n is a power of 2
		if(n != (1<<m))
			throw new RuntimeException("FFT length must be power of 2");

		if(cosHash.containsKey(n)){
			cos = (double[]) cosHash.get(n);
			sin = (double[]) sinHash.get(n);
		} else{
			// precompute tables
			cos = new double[n/2];
			sin = new double[n/2];

			for(int i=0; i<n/2; i++) {
				cos[i] = Math.cos(-2*Math.PI*i/n);
				sin[i] = Math.sin(-2*Math.PI*i/n);
			}
			cosHash.put(n,cos);
			sinHash.put(n,sin);
		}
	}

   /***************************************************************
  * fft.c
  * Douglas L. Jones 
  * University of Illinois at Urbana-Champaign 
  * January 19, 1992 
  * http://cnx.rice.edu/content/m12016/latest/
  * 
  *   fft: in-place radix-2 DIT DFT of a complex input 
  * 
  *   input: 
  * n: length of FFT: must be a power of two 
  * m: n = 2**m 
  *   input/output 
  * x: double array of length n with real part of data 
  * y: double array of length n with imag part of data 
  * 
  *   Permission to copy and use this program is granted 
  *   as long as this header is included. 
  ****************************************************************/  
	public void fft(double[] x, double[] y)
	{
		int i,j,k,n1,n2,a;
		double c,s,t1,t2;

		// Bit-reverse
		j = 0;
		n2 = n/2;
		for (i=1; i < n - 1; i++) {
			n1 = n2;
			while ( j >= n1 ) {
				j = j - n1;
				n1 = n1/2;
			}
			j = j + n1;

			if (i < j) {
				t1 = x[i];
				x[i] = x[j];
				x[j] = t1;
				t1 = y[i];
				y[i] = y[j];
				y[j] = t1;
			}
		}

		// FFT
		n1 = 0;
		n2 = 1;

		for (i=0; i < m; i++) {
			n1 = n2;
			n2 = n2 + n2;
			a = 0;

			for (j=0; j < n1; j++) {
				c = cos[a];
				s = sin[a];
				a +=  1 << (m-i-1);

				for (k=j; k < n; k=k+n2) {
					t1 = c*x[k+n1] - s*y[k+n1];
					t2 = s*x[k+n1] + c*y[k+n1];
					x[k+n1] = x[k] - t1;
					y[k+n1] = y[k] - t2;
					x[k] = x[k] + t1;
					y[k] = y[k] + t2;
				}
			}
		}
	}                          
}