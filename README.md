# multiSMD NAMD

This project is a user-friendly Python application that automates the generation of input files for multi-directional Steered Molecular Dynamics (SMD) simulations. The final application will be deployed as a web-based tool for streamlined execution.

Installation
--------------------------
This script, written in Python 3, relies on several external packages for its functionality. The required dependencies include NumPy, SciPy, and MDAnalysis.
The required dependencies can be installed in two ways:

1. Manual/implicit installation with pip or conda:

 ```pip install numpy scipy MDAnalysis```

```conda install numpy scipy MDAnalysis```

2. With provided `requirements.txt` file and _pip_ tool:

```pip install -r requirements.txt```

For consistent and reproducible results, the recommended method is to use the second one with provided requirements.txt file, which specifies the exact versions of all tested dependencies.

Runing
------------------------
Run the program in terminal (Linux) by typing in the command line (being in the directory where you downloaded the script):

_python multiSMD.py_   (or python3 multiSMD.py  if python3 is not the default python installation) in the same line we added parameters.

We call the program by giving it the following parameters:

_python multiSMD.py file.pdb file.psf file.vel file.coor file.xsc toppar.zip template.inp template.run ‘selection constraints’ ‘selection pull’_

file.pdb        - the PDB structure of our system

file.psf        - the PSF file for our system

file.vel        - velocity information- a restart file from NAMD equilibration simulation   

file.coor    - current coordinates of all atoms in the system (also from the NAMDs equilibration)

file.xsc        - contains the periodic cell parameters and extended system variables

toppar.zip    - if you have only one file with Charmm parameters, write it here
  (e.g. param1.inp). However, if there are more parameter files, you have to 
  put them into the toppar folder and "zip" it (_zip -r toppar.zip toppar_)

template.inp    - Here we have an input file to namd, in which we set all simulation parameters we need. Based on this file, the input to SMD will be generated, so __it is important that the section describing SMD and constraints (SMD on ... constraints yes, etc.) is present.__ The exemplary input file can be found in the TEST.zip folder. __To avoid possible artifacts, please make sure that SMD simulations are performed in the NVT ensemble__ (the pressure control should be set off).

template.run    - A sample script to run the simulation on the computer you intend to count - containing the namd running line. Input and output files will be defined as _INPF_ and _OUTF), so this is how they should be treated in the namd running line (_/home/user/NAMD/namd2 +p2 $INPF > $OUTF 2>&_). The exemplary template.run file can be found in TEST.zip.

‘selection’    - selections of constrained and pulled atoms 
 in SMD simulation. These are text variables, necessarily in quotation marks ''. The convention for atom selection is as in MDAnalysis (https://docs.mdanalysis.org/1.1.0/documentation_pages/selections.html) i.e., 'name CA and protein and segid A B C' or 'name CA and resid 1:55 66:128', or 'name CA and resname PRO ALA NBD'. Unfortunately, there is no 'chain' selection, so you have to use 'segid' instead.
It is recommended to restrain only CA atoms. The script produces an input file to Steered Molecular Dynamics (SMD_constraints.pdb) indicating which part should be restarined and which part will be pulled by changing values in columns O and B. 

The default mode of SMD simulation in here is constant Velocity SMD, which means that the pulling force will be adjusted to provide a constant velocity in the direction of pulling.

The output
--------------------------------
The program will generate an Output directory containing the input files and subdirectories corresponding to the "pull" directions. Each subdirectory contains an appropriately prepared input file to NAMD, and a bash script to run the given simulation (based on the provided template.run). You can run these scripts each separately:

 _. Output/SMD_theta_0_phi_0/run.bash_ (the "dot" will run the script exactly where the run.bash file is)
 
or together using the master.run script:

_./master.run_

These will start the SMD simulation. The default 

For visualisation purposes, the program generates a tcl script for VMD (_vmd_script.tcl_), which draw the bunch of vectors 
representing the directions of pulling. 

Analysing the output
---------------------------------
The Analysis.py script facilitates analyzing SMD simulation results. It iterates through subdirectories in the Output directory and extracts SMD information from mdrun.log files. Additionally, it checks the corresponding dcd files and calculates the number of hydrogen bonds formed between a given protein selection for each simulation.
The program can be called using the following parameters:

_python Analysis.py Output 'selection 1' 'selection 2' n_repeats_

where Output is the directory with SMD data, selection 1 is the constrained part of our system, and selection 2 is the pulled part of our system. 

# multiSMD GROMACS
To run multiSMD simulations in gromacs you also need to start from pre-equilibrated system input.
Run the program in terminal (Linux) by typing in the command line (being in the directory where you downloaded the script):

_python multiSMD_GRO.py_   (or python3 multiSMD_GRO.py  if python3 is not the default python installation) in the same line we added parameters.

We call the program by giving it the following parameters:

_python multiSMD.py file.pdb file.gro restart_md.zip template.mdp ‘selection constraints’ ‘selection pull’ n_repeats_

file.pdb        - the PDB structure of our system (it is needed here for proper chain name information)

file.gro        - the GRO file for our system

restart_mdr.zip    - all files necessary to restart your simulation - topol.top, toppar, index.ndx etc. (including all forcefield parameters)

template.mdp    - The template mdp file needs to have the SMD-controlling fragment already inside The exemplary input file can be found in the TEST_GRO.zip folder. __To avoid possible artifacts, please make sure that SMD simulations are performed in the NVT ensemble__ (the pressure control should be set off).

‘selection’    - selections of constrained and pulled atoms 
 in SMD simulation. These are text variables, necessarily in quotation marks ''. The convention for atom selection is as in MDAnalysis (https://docs.mdanalysis.org/1.1.0/documentation_pages/selections.html) i.e., 'name CA and protein and segid A B C' or 'name CA and resid 1:55 66:128', or 'name CA and resname PRO ALA NBD'. Unfortunately, there is no 'chain' selection, so you have to use 'segid' instead.
It is recommended to restrain only CA atoms. The script produces an input file to Steered Molecular Dynamics (SMD_constraints.pdb) indicating which part should be restarined and which part will be pulled by changing values in columns O and B. 

n_repeats - how many repeats of each pulling direction you wan to perform

The default mode of SMD simulation in here is constant Velocity SMD, which means that the pulling force will be adjusted to provide a constant velocity in the direction of pulling.

The program will generate an Output directory containing the input files and subdirectories corresponding to the "pull" directions. Each subdirectory contains an appropriately prepared input file to GROMACS.

Analysing the output
---------------------------------
The Analysis_GRO.py script facilitates analyzing SMD simulation results. It iterates through subdirectories in the Output directory and extracts SMD information from mdrun.log files. Additionally, it checks the corresponding xtc files and calculates the number of hydrogen bonds formed between a given protein selection for each simulation.
The program can be called using the following parameters:

_python Analysis_GRO.py Output 'selection 1' 'selection 2' n_repeats file.pdb_

where Output is the directory with SMD data, selection 1 is the constrained part of our system, and selection 2 is the pulled part of our system. 
