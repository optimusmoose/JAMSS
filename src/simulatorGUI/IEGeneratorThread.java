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

import java.util.ArrayList;
import java.util.LinkedList;

/**
 *
 * @author rob
 */
public class IEGeneratorThread extends Thread{
	ArrayList<Peptide> queue;
	private int threadId;
	public boolean finished = false;
	public MassSpec massSpec;
	public IEGeneratorThread(ArrayList<Peptide> _queue, int threadId){
		queue = _queue;
		massSpec = new MassSpec(threadId);
	}
	
	@Override
	public void run() {
		// create one mass spec object per thread
		for(Peptide peptide:queue) {
      // write out to disk once memory is scant
			if(Runtime.getRuntime().freeMemory() < Runtime.getRuntime().maxMemory()*0.10){
         serializeEnvelopes();
			}
			//if(peptide==null){return;}
			processPeptide(peptide);
      
      LocalProgressMonitor.updateIsotopePatternSimulation();
		}
		serializeEnvelopes();
    massSpec.isotopicEnvelopesInstance = null;
		finished = true;
		return;
	}
  
  private void serializeEnvelopes(){
    for(IsotopicEnvelope ie : massSpec.isotopicEnvelopesInstance){
      ie.serialize(MassSpec.pathToClass,threadId);
    }
    massSpec.isotopicEnvelopesInstance = new LinkedList<IsotopicEnvelope>();
  }
  
  public boolean processPeptide(Peptide peptide){
		if (MassSpec.modifications.powerSet.size() == 0){ // no variable mods
			massSpec.computeIsotopicEnvelopes(peptide, new ArrayList<Modification>());
		} else {
      double totalAbundance = peptide.peptideIntensity;
      for (int i=0; i < MassSpec.modifications.powerSet.size(); i++){ // for each combination of mods
				// create peptide with the proper intensity
        peptide.peptideIntensity = totalAbundance * ((Modification) MassSpec.modifications.powerSet.get(i).get(0)).percent;
				massSpec.computeIsotopicEnvelopes(peptide, MassSpec.modifications.powerSet.get(i));
			}
    }
		return true;
	}
}


