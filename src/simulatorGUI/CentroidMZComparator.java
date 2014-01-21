/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package simulatorGUI;

import java.util.Comparator;

/**
 *
 * @author rob
 */
public class CentroidMZComparator implements Comparator<Centroid>  {
    @Override
    public int compare(Centroid o1, Centroid o2) {
        return ((Double) o1.getMZ()).compareTo((Double) o2.getMZ());
    }
}
