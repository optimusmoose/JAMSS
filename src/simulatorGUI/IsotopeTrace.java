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
	private static final double TWOPI = Math.PI * 2.0;
	private static final double EXPZERO = Math.exp(0.0000000000000001);
	
	double sdShareY;
	double sdShareZ;
	double sdX;
	double sdY;
	double sdZ;
	double centroidIntensity;
	double oneDIntensityFactor;
	double maxXIntensity;
	double normalizingConstantY;
	double normalizingConstantZ;
	double muGaussian;
	double sdShareX;
	double sdEstimate;
	double predictedIEIntensity;
	double isotopeTraceIntensityRatio;
	double traceLength;
	double predictedRt;
	double isotopeTraceMass;
	public IsotopeTrace(double _sdShareX, 
						double _sdEstimate, 
						double _predictedIEIntensity, 
						double _isotopeTraceIntensityRatio, 
						double _traceLength, 
						double _predictedRt, 
						double _isotopeTraceMass){
		sdShareX = _sdShareX;
		sdEstimate = _sdEstimate;
		predictedIEIntensity = _predictedIEIntensity;
		isotopeTraceMass = _isotopeTraceMass;
		isotopeTraceIntensityRatio = _isotopeTraceIntensityRatio;
		traceLength = _traceLength;
		predictedRt = _predictedRt;
		sdShareY = 0.35 + RandomFactory.rand.nextDouble() * 0.15; // random between 0.35-0.50
		sdShareZ = 1.0 - (sdShareY + sdShareX);
		sdX = sdShareX * sdEstimate;
		sdY = sdShareY * sdEstimate;
		sdZ = sdShareZ * sdEstimate;

		centroidIntensity=0;
		maxXIntensity = predictedIEIntensity * isotopeTraceIntensityRatio * (1.0/(sdX*Math.sqrt(TWOPI))) * EXPZERO;
		normalizingConstantY = maxXIntensity / (2*(1.0/(sdY*Math.sqrt(TWOPI))) * EXPZERO);
		normalizingConstantZ = maxXIntensity / (2*(1.0/(sdZ*Math.sqrt(TWOPI))) * EXPZERO);
		muGaussian = predictedRt + traceLength*sdShareX;

	}
	
	public Centroid getCentroidAtRT(double rt){
		Centroid result = new Centroid(0,0);
		if (MassSpec.oneD) {
			oneDIntensityFactor = 0.05 + 0.45 * RandomFactory.rand.nextFloat(); // between 0.05 and 1
			centroidIntensity = predictedIEIntensity * isotopeTraceIntensityRatio * oneDIntensityFactor;
		} else { // normal chromotography (not 1 d)
			double gaussianX = predictedIEIntensity * isotopeTraceIntensityRatio * (1.0/(sdX*Math.sqrt(TWOPI))) * Math.exp(-(Math.pow(rt-muGaussian,2)/(2.0 * Math.pow(sdX,2))));
			double gaussianY = normalizingConstantY * (1.0/(sdY*Math.sqrt(TWOPI))) * Math.exp(-(Math.pow(rt-muGaussian,2)/(2.0 * Math.pow(sdY,2)))); 
			double gaussianZ = normalizingConstantZ * (1.0/(sdZ*Math.sqrt(TWOPI))) * Math.exp(-(Math.pow(rt-muGaussian,2)/(2.0 * Math.pow(sdZ,2)))); 
			if (rt < predictedRt + sdShareX * traceLength){ // first half treated as a more narrow Gaussian
				centroidIntensity = gaussianX + (gaussianY + gaussianZ)*(gaussianX / maxXIntensity);
			} else { // second half treated as a wider Gaussian
				centroidIntensity = gaussianX + gaussianY + gaussianZ;
			}
		}
		if (centroidIntensity > MassSpec.minWhiteNoiseIntensity){ //filter points that are too small
			// jaggedness - insert predictedIntensity noise
			centroidIntensity -= (centroidIntensity/2.0) * RandomFactory.rand.nextDouble() * RandomFactory.rand.nextDouble();
			// mz wobble - insert mz noise. Two types: Universal and intensity specific. 
			// universal:
			double universalNoised = isotopeTraceMass - 0.001 + RandomFactory.rand.nextDouble() * 0.002; //rand amt between -0.001 and 0.001
			// intensity-specific:
			double maxWobble = 0.03 * (MassSpec.minWhiteNoiseIntensity / centroidIntensity);
			// get random noise amount between min mz and max mz
			double mzWobble = universalNoised + (RandomFactory.rand.nextBoolean() ? -1 : 1)* (RandomFactory.rand.nextDouble() * maxWobble);
			if (mzWobble > 0){
				// TADA! finished centroid. Add to output map.		
				result = new Centroid(mzWobble, centroidIntensity);
				if (mzWobble > MassSpec.maxMZ){MassSpec.maxMZ = mzWobble;}
			}
		}
		return result;
	}
}
