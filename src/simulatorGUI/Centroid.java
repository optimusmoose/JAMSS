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
public class Centroid implements java.io.Serializable{
	public double mz;
	public double abundance;
	public static int totalCentroids=0;
	public int centroidID;
	public int ionFeatureID;
	public int charge;
	public int pepID;
	public int proteinID;
	
	private Centroid(){}
	
	public Centroid(double _mz, double _abundance){
		centroidID = totalCentroids;
//System.out.println(totalCentroids);
		totalCentroids += 1;
		mz = _mz;
		abundance = _abundance;
	}
	
	public double getMZ(){
		return mz;
	}
	
	public Centroid clone(){
		Centroid cent = new Centroid();
		cent.abundance = this.abundance;
		cent.centroidID = this.centroidID;
		cent.charge = this.charge;
		cent.ionFeatureID = this.ionFeatureID;
		cent.mz = this.mz;
		cent.pepID = this.pepID;
		cent.proteinID = this.pepID;
		return cent;
	}
}
