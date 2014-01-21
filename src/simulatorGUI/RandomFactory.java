/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package simulatorGUI;

import java.util.Random;

/**
 *
 * @author rob
 */
public class RandomFactory {
	// TODO fold in digester options
	public Random rand;
	public long randSeed;
	public RandomFactory() {
		randSeed = System.currentTimeMillis();
		rand = new Random();
		rand.setSeed(randSeed);
	}
	public RandomFactory(long _randseed){
		randSeed = _randseed;
		rand = new Random();
		rand.setSeed(randSeed);
	}
}
