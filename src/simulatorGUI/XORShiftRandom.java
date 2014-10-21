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

import java.util.Random;

/**
 *
 * BASED ON @article{marsaglia2003xorshift,
  title={Xorshift rngs},
  author={Marsaglia, George},
  journal={Journal of Statistical Software},
  volume={8},
  number={14},
  pages={1--6},
  year={2003}
}
 */
public class XORShiftRandom extends Random{
	public long seed;

	public XORShiftRandom() {
		seed = System.nanoTime();
	}
	@Override public void setSeed(long _seed){
		seed = _seed;
	}
	@Override protected int next(int nbits) {
		// N.B. Not thread-safe!
		long x = seed;
		x ^= (x << 21);
		x ^= (x >>> 35);
		x ^= (x << 4);
		seed = x;
		x &= ((1L << nbits) -1);
		return (int) x;
	}
}
