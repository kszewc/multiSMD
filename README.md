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

Launching SMD scripts
------------------------
Run the SMD scripts in terminal by typing in the command line:

- For NAMD: ``` multismd_namd.py [-h] [--repeats REPEATS] output_dir input_pdb input_psf input_vel input_coor input_xsc input_par1 template_inp template_run input_sel1 input_sel2```   

  where:
    - output_dir:         The main output directory for all generated files.
    - input_pdb:          Path to the NAMD MD restart PDB file. 
    - input_psf:          Path to the NAMD MD restart PSF file.
    - input_vel:          Path to the NAMD MD restart velocity file.
    - input_coor:         Path to the NAMD MD restart coordinate file.
    - input_xsc:          Path to the NAMD MD restart extended system file.
    - input_par1:         Path to the primary force field parameter file.
    - template_inp:       Path to the template NAMD input file (.inp).
    - template_run:       Path to the template bash run script (.run).
    - input_sel1:         MDAnalysis selection criteria for the constrained atoms (fixed group).
    - input_sel2:         MDAnalysis selection criteria for the pulled atoms (pulled group).

- For Gromacs: ```multismd_gromacs.py [-h] [--repeats REPEATS] output_dir input_pdb input_gro input_md template_mdp input_sel1 input_sel2```
 
  where:
    - output_dir:         The main output directory for all generated files.
    - input_pdb:          Path to the MD PDB file.
    - input_gro:          Path to the Gromacs GRO file.
    - input_md:           Zipped Gromacs files needed to restart simulation (includes: .gro .top toppar .ndx)
    - template_mdp:       Paths to Gromacs input file and run.bash file (.mdp)
    - input_sel1:         MDAnalysis selection criteria for the constrained atoms (fixed group).
    - input_sel2:         MDAnalysis selection criteria for the pulled atoms (pulled group).

**Detailed parameters descriptions** :


**input_par1**    - if you have only one file with Charmm parameters, write it here
  (e.g. param1.inp). However, if there are more parameter files, you have to 
  put them into the toppar folder and "zip" it (_zip -r toppar.zip toppar_)

**template.inp**    - Based on this file, the input to SMD will be generated, so __it is important that the section describing SMD and constraints (SMD on ... constraints yes, etc.) is present.__ The exemplary input file can be found in the TEST.zip folder. __To avoid possible artifacts, please make sure that SMD simulations are performed in the NVT ensemble__ (the pressure control should be set off).

**template.run**    - A sample script to run the simulation on the computer you intend to count - containing the namd running line. Input and output files will be defined as _INPF_ and _OUTF), so this is how they should be treated in the namd running line (_/home/user/NAMD/namd2 +p2 $INPF > $OUTF 2>&_). The exemplary template.run file can be found in TEST.zip.

**input_sel1** and **input_sel2**    - selections of constrained and pulled atoms 
 in SMD simulation. These are text variables, necessarily in quotation marks ''. The convention for atom selection is as in MDAnalysis (https://docs.mdanalysis.org/1.1.0/documentation_pages/selections.html) i.e., 'name CA and protein and segid A B C' or 'name CA and resid 1:55 66:128', or 'name CA and resname PRO ALA NBD'. Unfortunately, there is no 'chain' selection, so you have to use 'segid' instead.
It is recommended to restrain only CA atoms. The script produces an input file to Steered Molecular Dynamics (SMD_constraints.pdb) indicating which part should be restarined and which part will be pulled by changing values in columns O and B. 

The default mode of SMD simulation in here is constant Velocity SMD, which means that the pulling force will be adjusted to provide a constant velocity in the direction of pulling.

Output
--------------------------------
The program generates by a directory containing the input files and subdirectories corresponding to the "pull" directions. 
Each subdirectory contains an appropriately prepared input file to NAMD or Gromacs, and a bash script to run the given simulation (based on the provided template.run). 
For NAMD you can run these scripts each separately:

```. Output/SMD_theta_0_phi_0/run.bash``` (the "dot" will run the script exactly where the run.bash file is)
 
or together using the master.run script:

```./master.run```

These will start the SMD simulation.

For visualisation purposes, the program generates a tcl script for VMD (_vmd_script.tcl_), which draw the bunch of vectors 
representing the directions of pulling. 

Analysing the output
---------------------------------
The Analysis.py script facilitates analyzing SMD simulation results. 
It iterates through subdirectories in the output directory and extracts SMD information from mdrun.log files. Additionally, it checks the corresponding dcd files and calculates the number of hydrogen bonds formed between a given protein selection for each simulation.
The program can be called using the following parameters:

- For NAMD:
```python analysis_namd.py output_dir 'selection 1' 'selection 2' n_repeats```

- For Gromacs:
```python analysis_gromacs.py output_dir 'selection 1' 'selection 2' n_repeats```

where:
 - output_dir - is the directory with SMD data
 - selection 1 - is the constrained part of our system
 - selection 2 - is the pulled part of our system