/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package simulatorGUI;

import java.util.ArrayList;

/**
 *
 * @author rob
 */
public class TestDigester {
	public static void main(String [ ] args){
		Digester digester = DigesterOpts.getDigester("trypsin");
		digester.testProcessProtein();
		ArrayList<String> result = digester.processProtein("HLKTEAEMK", 2);
		for(String s : result){
			System.out.println(s);
		}
	}
}
