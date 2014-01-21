		JAMSS is a mass spec simulator.

		Copyright (C) 2014  Rob Smith (2robsmith@gmail.com)

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

Minimumalist run:
1) Open a fasta file. (IMPORTANT: Make sure fasta headers do not contain a #. If they do, remove it.)
2) Click run.

How to clone a run: 
1) Open an mzML file created by JaMMS with the "Load Options" Open button. All parameters will be automatically populated in the GUI.
2) Unclick the "Replicate?" box. This means you want to use the same random seed as was used in the original file, and you will get the same exact output file as you loaded options from, provided the .fasta file has not changed.
3) Click run.

How to replicate a run:
1) Open an mzML file created by JaMMS with the "Load Options" Open button. All parameters will be automatically populated in the GUI.
2) Optional: tweak options. The more things you change, the more different the "replicate" will be.
3) Hit run.

More involved run:
1) Open a fasta file.
2) Turn mods on/off by selecting check box (static modifications) or adjusting percentage affected (variable modifications). NOTE: using variable modifications will significantly increase runtime, as it increases the number of simulation permutations that have to be calculated.
3) Adjust how many cores you want to use. If you want to continue working on this machine at the same time, select at least one less than the total number available. 
4) Adjust the digestion enzyme.
5) Adjust the scans per second. The higher the number, the longer the program will take to run.
6) Adjust the runtime. The higher the number, the longer the program will take to run.
7) Adjust the missed cleavages. They will include all permutations less than the number selected, so the size of the simulation increases substantially with each increase.
8) Select the dropout percentage. This is how many scans will randomly be excluded from the simulation.
9) Adjust the white noise points.
10) The resolution of the simulation can be adjusted by selecting a different overlap range. This value is in mz.
11) The intensities of each protein can be set in the fasta file. For each header, append a decimal value indicating the percent of the sample composed of this protein preceded by a pound sign. For example, #0.342. The percentages of the whole file should sum to 1.
12) Hit run.

How to simulate a non-chromatographic run:
1) Follow above instructions, but check the "One Dimension Simulation" box.
2) Hit run.

Advanced options to take caution with:
1) Overlap Range. This controls the effective resolution of the mass spectrometer (higher number, lower resolution).


Ways to make faster runs:
1) Decrease modifications (significant).
2) Increase number of CPUs (significant).
3) Increase amount of memory allocated to Java Virtual Machine. This can be done by adding "Xmx10g" onto the command line of your java run, where 10 is the number of GB of RAM you are allocating to the JVM (significant).
4) Decrease the number of missed cleavages (significant).
5) Decrease the number of MS2 per scan (minor).
6) Decrease number of scans per second, noise points, or runtime (minor).
7) Increase minimum noise intensity, which is the global floor of points that are accepted in simulation (minor).
8) Increase overlap range (careful with this one) (minor).




