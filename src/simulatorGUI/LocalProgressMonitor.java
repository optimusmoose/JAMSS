/*
 * Copyright (C) 2014 rob
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
 * Class to organize progress bar. Starts anticipating time for digestion.
 * @author rob
 */
public class LocalProgressMonitor {
  static int totalComputationEpochs;
  static int currComputationEpochs;
  static int currIdx;
  static private long lastTime;
  
  public static void checkCancel(){
    if (simulatorGUI.progressMonitor.isCanceled()) {
      System.exit(1);
    }
  }
  
  public static void beginDigestion(int numFastaLines){
    checkCancel();
    lastTime = System.nanoTime();
    totalComputationEpochs = numFastaLines;
    simulatorGUI.progressMonitor.setNote("Phase 1/3: Digesting...");
		simulatorGUI.progressMonitor.setProgress(0);
  }
  
  public static void updateDigestion(int fastaIndex){
    checkCancel();
    long now = System.nanoTime();
    long timePerFasta = now - lastTime;
    lastTime = now;
    simulatorGUI.progressMonitor.setNote("Phase 1/3: Digesting...(" + 
    (int) ((totalComputationEpochs-fastaIndex) * timePerFasta * 0.000000000016667) + " mins remaining");
		simulatorGUI.progressMonitor.setProgress((int) (100.0 * ((double) fastaIndex / (double) totalComputationEpochs)));
  }
  
  public static void beginIsotopePatternSimulation(int numPeptides){
    checkCancel();
    lastTime = System.nanoTime();
    totalComputationEpochs = numPeptides;
    simulatorGUI.progressMonitor.setNote("Phase 2/3: Simulating Isotope Patterns...");
		simulatorGUI.progressMonitor.setProgress(0);
    currIdx = 0;
  }
  
  public static void updateIsotopePatternSimulation(){
    checkCancel();
    long now = System.nanoTime();
    long timePerScan = now - lastTime;
    lastTime = now;
    if(currIdx % 10 == 0){
      simulatorGUI.progressMonitor.setNote("Phase 2/3: Simulating Isotope Patterns..." +
        (int) ((totalComputationEpochs-currIdx) * timePerScan * 0.000000000016667) + " mins remaining");
      simulatorGUI.progressMonitor.setProgress((int)(100.0 * ((double) currIdx / (double) totalComputationEpochs)));
    }
    currIdx++;
  }
  
  public static void beginExtendSimulation(int numScans){
    checkCancel();
    currIdx = 0;
    lastTime = System.nanoTime();
    totalComputationEpochs = numScans;
    simulatorGUI.progressMonitor.setNote("Phase 3/3: Extending Isotope Patterns...");
		simulatorGUI.progressMonitor.setProgress(0);
  }
  
  public static void updateExtendSimulation(int scanIndex){
    checkCancel();
    long now = System.nanoTime();
    long timePerScan = now - lastTime;
    lastTime = now;
    simulatorGUI.progressMonitor.setNote("Phase 3/3: Extending Isotope Patterns..." +
      (int) ((totalComputationEpochs-scanIndex) * timePerScan * 0.000000000016667) + " mins remaining");
		simulatorGUI.progressMonitor.setProgress((int)(100.0 * ((double) scanIndex / (double) totalComputationEpochs)));
    currIdx++;
  }
  
}
