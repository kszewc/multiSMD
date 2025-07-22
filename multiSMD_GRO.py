###########################################################################################
#
#   Script creating an input for the multi-directional SMD simulation in GROMACS. 
#   Series of Steered MD simulations are performed each in different direction of pulling 
#   to roughly cover a semisphere. Each simulation data is stored in separate directory. 
#   To run it, you need to have numpy, scipy and MDAnalysis libraries installed.
#
#   Author: Katarzyna Walczewska-Szewc, Nicolaus Copernicus University in Toruń, 06.12.2024
###########################################################################################

from MDAnalysis import *
from numpy import *
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
import os
import sys
from zipfile import ZipFile
##########################################################################################
#
#   The input data and parameters
#
##########################################################################################

dirname = 'Output' # name of the new directory
# GROMACS MD restarting files
input_pdb = str(sys.argv[1]) #'it is important to use pdb to have information about chain names'
input_gro = str(sys.argv[2]) #actual frame to start with
input_md = str(sys.argv[3]) #files necesaery to restart md simulations in gromacs: .gro .top toppar .ndx - all zipped, includes forcefield parameters

# Templates for gromacs input file and run.bash file
template_mdp = str(sys.argv[4]) #'../CUT/template.mdp'
#template_run = str(sys.argv[3]) #'../CUT/template.run' #not implemented for GROMACS - will appear in the future

#Selection parameters for SMD (see MDAnalysis page for more information about selection keywords)
input_sel1= sys.argv[5] # atoms to pull something from (co calculate proper direction)
input_sel2= sys.argv[6] # atoms to be pulled

n_repeats = 1# int(sys.argv[11]) # number of simulation repeats for each pulling direction

###########################################################################################
#
#   End of parameters setup, do not change unless you know what you are doing ;)
#
###########################################################################################


# Function generating the set of vectors samplig hemisphere. 
# Such hemisphere is alighed to the direction of the vector joining COM of the fix and pull selections from the pdb file.

def Hedgehog(pdb,actual,fix_selection,pull_selection):

    u = Universe(pdb,actual)

    fix = u.select_atoms(fix_selection)
    pull = u.select_atoms(pull_selection)
    fix_COM = fix.center_of_mass()
    pull_COM = pull.center_of_mass()

    ax_principal = pull_COM - fix_COM #creating vector pointing the general direction of the cone
    ax_principal = ax_principal/linalg.norm(ax_principal)#creating the unit vector
    print((ax_principal))

    #Generating a hedgehog of vectors in z-direction
    hedgehog = array([[0,0,1]])
    labels = array([[0, 0]])
    for theta,resolution in zip([45,90],[90,90]):#[15,30,45,60,75,90],[90,75,60,45,30,15]):
        for phi in range(0,360,resolution):
            x = cos(deg2rad(phi))*sin(deg2rad(theta))
            y = sin(deg2rad(phi))*sin(deg2rad(theta))
            z = cos(deg2rad(theta))

            hedgehog = concatenate((hedgehog,array([[x,y,z]])),axis=0)
            labels = concatenate((labels,array([[theta,phi]])),axis=0) #labels for subdirectories names construction

    #Transforming the vectors cone to the direction given by the principal ax:
    r = R.align_vectors([ax_principal],[[0,0,1]])
    hedgehog = r[0].apply(hedgehog)
    

    return hedgehog, labels, pull_COM #set of vectors [N][x,y,z]; set of angles [N][theta,phi], the pulling point [x,y,z]

##############################################################################################
#
#   Function generating a tree of directories containing input and run  files for each SMD run
#   The template.inp and run.bash files are modified and placed in proper directories. The copies of the GROMACS restart files
#   are created in the output directory, to make it self-sustainable.
#
##############################################################################################




def Gen_input(name,pdb,vectors,template,sel1,sel2,mdf):
    
    #   name/ - the output directory containing all generated files, copy of the NAMD input structures and the simulation results
    #       SMD_constraints.pdb - the output file containing information about fixed and pulled atoms for SMD (O and B colum, respectively)
    #       SMD_theta_i_phi_j/ - subdirectories for each SMD direction run (the simulation output will be stored here)
    #           mdrun.mdp - input file for a single SMD simulation
    #           run.bash - the run.bash script (generated in next function)
    #       vmd_script.tcl - the vmd script to visualize the cone of vectors (generated below)
    
    # Creating copy of GROMACS restart files in the output directory
   
    os.system(f'rm -rf {name}')
    os.system('mkdir '+str(name))
    os.system('cp '+str(pdb)+' '+str(name)+'/')


    # Creating SMD_constraints.pdb file with flags in O and B column indicating fixed and pulled atoms
    u = Universe(pdb)

    pul = u.select_atoms(sel2)


    # Creating subdirectories for all directions in the semisphere
    for v,l in zip(vectors[0],vectors[1]):
        print('theta_' + str(int(l[0])) + '_phi_' + str(int(l[1])))

        with ZipFile(mdf,'r') as z:
             z.extractall(str(name)+'/')
        m = mdf.split('/')[-1]
        os.system(f'mv {name}/{m[:-4]} ' + str(name) + '/SMD_theta_' + str(int(l[0])) + '_phi_' + str(int(l[1])))

        # Modyfying the mdrun.mdp file for each SMD run (the only difference between them is the direction of pulling)
        

        f = open(template,'r')
        new_f= open(f'{name}/SMD_theta_{int(l[0])}_phi_{int(l[1])}/mdrun.mdp','w')  # several repearts of simulation is produced


        for line in f.readlines():
            if line[0:15] == 'pull-coord1-vec':
                new_f.write(f'pull-coord1-vec		= {v[0]} {v[1]} {v[2]} \n')
            elif line[0:15] == 'pull-group1-pbc':
                new_f.write(f'pull-group1-pbcatom	= {pul.atoms[0].index} \n')
            else:
                new_f.write(line)
        new_f.close()
        f.close()

##############################################################################################
#
#   Function generating run.bash files inside the tree of subdirectories for each SMD run
#   This function also initiate the simulation so if you do not want it to run immediately,
#   please comment out the proper line. The template.bash file should be modified to 
#   fit your computing facility.
#
##############################################################################################



def Gen_run(name,label,template): #not working yet for gromacs - please run your simulations manually
    master_file = open(str(name)+'/master.run','w')
    for l in label:
        # Generating run.bash files
        
        for n in range(0,n_repeats):
                f = open(template,'r')
                new_f= open(str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+f'/run{n}.bash','w')
                for line in f.readlines():
                    line_new = line
                    if len(line.split())>1 and line.split()[1] == '-J': line_new = line.split()[0] + ' ' + line.split()[1] +' SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+'\n'
                    
                    if len(line.split())>1 and line.split()[1][0:4] == 'INPF': line_new = 'set INPF=SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+f'/mdrun{n}.inp'+'\n'
                    if line[0:4] == 'INPF': line_new = f'INPF=mdrun{n}.inp'+'\n'
                    if line[0:4] == 'OUTF': line_new = f'OUTF=mdrun{n}.log'+'\n'

                    if len(line.split())>1 and line.split()[1][0:4] == 'OUTF': line_new = 'set OUTF=SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+f'/mdrun{n}.log'+'\n'
                    new_f.write(line_new)
                #print('qsub '+str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+f'/run{n}.bash')
                new_f.close()
        

                # Runing the simulations
                os.system('chmod +x '+str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+f'/run{n}.bash')
                #print('Running simulation for SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1])))
                #os.system('. '+str(name)+'/SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+f'/run{n}.bash &') # If you are using this script on a supercomputer, all jobs can run simultanously 
                master_file.write('. SMD_theta_'+str(int(l[0]))+'_phi_'+str(int(l[1]))+f'/run{n}.bash\n') # On a single gpu station, it could be better to run one job after one (the line without & at the end so the next job is waiting for the one to finish)
                f.close()
    master_file.close()
    os.system('chmod +x '+str(name)+'/master.run')
        
#########################################################################################
#
#   The tcl script that can be run in the VMD tcl console to visualize the cone of vectors
#
#########################################################################################
def Gen_vmd_script(name,vector,COM):
    f = open(name+'/vmd_script.tcl','w')
    #array draving subroutine (from https://www.ks.uiuc.edu/Research/vmd/current/ug/node127.html)
    f.write('proc vmd_draw_arrow {mol start end color label} {\n\t# an arrow is made of a cylinder and a cone\n\t set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]\n\t graphics $mol color $color\n\t graphics $mol cylinder $start $middle radius 0.15\n\t graphics $mol cone $middle $end radius 0.25\n')
    f.write('set text_pos [vecadd $middle [list 0.5 0.5 0.5]] ; # offset text slightly for visibility\n')
    f.write('graphics $mol text $text_pos $label size 1\n}\n')
    s=5 #scaling the vectors
    nr_vec=0
    for j in vector:
        i = j + COM # The vectors need to be translated to the center of the mass of the protein and rescaled by the scaling factor (below)
        c = 1 #vmd color 1=blue
        f.write('draw arrow {'+str(COM[0])+' '+str(COM[1])+' '+str(COM[2])+'} {'+str(i[0]+(i[0]-COM[0])*s)+' '+str(i[1]+(i[1]-COM[1])*s)+' '+str(i[2]+(i[2]-COM[2])*s)+'} 1 \"'+str(nr_vec)+'\"\n' )
        nr_vec += 1




#########################################################################################
#
#   The body of the script
#
#########################################################################################

hedgehog = Hedgehog(input_pdb,input_gro, input_sel1, input_sel2)
Gen_input(dirname, input_pdb, hedgehog, template_mdp, input_sel1, input_sel2, input_md)
#Gen_run(dirname, hedgehog[1], template_run) #not working yet for GROMACS
Gen_vmd_script(dirname, hedgehog[0], hedgehog[2])
print('Generation of SMD input finished!')
