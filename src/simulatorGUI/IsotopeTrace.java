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

/**
 *
 * @author Rob.Smith
 */
public class IsotopeTrace {
	int rtIndex;
	double centroidIntensity;
	double oneDIntensityFactor;
	double maxXIntensity;
	double normalizingConstantY;
	double normalizingConstantZ;
	int traceLength; // in terms of number of scans
	double predictedRt;
	double isotopeTraceMass;
	double rand1;
	public IsotopeTrace(double _predictedIEIntensity, 
						double _isotopeTraceIntensityRatio, 
						int _traceLength, 
						double _predictedRt, 
						double _isotopeTraceMass,
						double _rand1){
		isotopeTraceMass = _isotopeTraceMass;
		traceLength = _traceLength;
		predictedRt = _predictedRt;
		rtIndex = 1;
		rand1 = _rand1 * 0.20;
		centroidIntensity=0;
		maxXIntensity = _predictedIEIntensity * _isotopeTraceIntensityRatio;

	}
	
	public Centroid getCentroidAtRT(double rt, RandomFactory localRandomFactory){
		Centroid result = new Centroid(0,0);
		if (MassSpec.oneD) {
			oneDIntensityFactor = 0.05 + 0.45 * localRandomFactory.localRand.nextFloat(); // between 0.05 and 1
			centroidIntensity = maxXIntensity * oneDIntensityFactor;
		} else { // normal chromotography (not 1 d)
			//centroidIntensity = maxXIntensity * Math.exp(-Math.pow(rtIndex - predictedRt,2.0) / ( ((double) rtIndex / (double) traceLength) * (traceLength + 14.66692)));
			centroidIntensity = maxXIntensity * Math.exp(-Math.pow(rtIndex - predictedRt,2.0) / Math.pow(((double) rtIndex / (double) traceLength) * traceLength * 0.14,2.0));
		}

		if (centroidIntensity > MassSpec.minWhiteNoiseIntensity){
			double wobbleJaggedRand = localRandomFactory.localRand.nextDouble();
			double minWobble = 0.001;
			double maxWobble = 0.02;
			double maxIntScaled = 5.0 / MassSpec.maxIntensity * centroidIntensity;
			double mzWobble = Math.max(minWobble,maxWobble * Math.exp(-maxIntScaled));
			mzWobble = isotopeTraceMass - mzWobble + mzWobble * 2.0 * wobbleJaggedRand; 
				// TADA! finished centroid. Add to output map.		
				result = new Centroid(mzWobble, centroidIntensity - 0.30 * centroidIntensity + 0.60 * centroidIntensity * wobbleJaggedRand);
		}
		rtIndex++;
		return result;
	}
}
