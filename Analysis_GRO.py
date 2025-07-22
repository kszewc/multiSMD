from numpy import *
from MDAnalysis import *
import os,sys
import matplotlib.pyplot as plt
import matplotlib
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from mpl_toolkits.mplot3d import Axes3D



def Dist(a, b):
    r = sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))
    return r

def plot_vectors(filename, theta,phi, ax, start_point=(0, 0, 0)):
    # Draw a bunch of vectors from file and highlight one to red


    # loading data
    data = loadtxt(filename)

    # Setting a vector list
    vectors = data - start_point

    # highlighted vector
    x = cos(deg2rad(phi))*sin(deg2rad(theta))
    y = sin(deg2rad(phi))*sin(deg2rad(theta))
    z = cos(deg2rad(theta))
    red = [x,y,z]

    # Drawing vectors
    for vector in vectors:
        ax.quiver(start_point[0], start_point[1], start_point[2],
                  vector[0], vector[1], vector[2],
                  length=1.5, normalize=True, color='gray')
    ax.quiver(start_point[0], start_point[1], start_point[2],
                  red[0], red[1], red[2],
                  length=1.5, normalize=True, color='red')

    #ax.text(start_point[0] + red[0]*1.1, start_point[1] + red[1]*1.1, start_point[2] + red[2]*1.1,
    #                f"θ={int(theta)}, φ={int(phi)}", color='red')
    ax.set_title(f"θ={int(theta)}, φ={int(phi)}")
    # Plot range
    ax.set_xlim(min(data[:, 0]) - 0.1, max(data[:, 0]) + 0.1)
    ax.set_ylim(min(data[:, 1]) - 0.1, max(data[:, 1]) + 0.1)
    ax.set_zlim(min(data[:, 2]) - 0.1, max(data[:, 2]) + 0.1)
    ax.set_axis_off()

name = str(sys.argv[1]) # Name of your SMD output directory
pdb = str(sys.argv[5]) #pdb file for structure related calculation
sel_const = str(sys.argv[2]) # First selection (to calculate center of mass and H-bonds) - constrained part
sel_pull = str(sys.argv[3]) # Second selection (to calculate center of mass and H-bonds) - pulled part
n_repeats = int(sys.argv[4]) # Here is the number of simulation repeats for each pulling direction
#extract_forces = False # change it to False after force files are generated
calculate_HB = False # change it to False after HB files are generated
red = 5 # reducing factor for too dense data

list_of_dirs = [subdir for subdir in os.listdir(name) if subdir[0:9] == 'SMD_theta']

## this part takes some time, since it is analysing your SMD trajectories looking for possible H-bonds between selections
u = Universe(name+'/'+pdb)
fix = u.select_atoms(sel_const).center_of_mass() #
bunch = [] #here the information about pulling vectors will be stored

for dir in list_of_dirs:
    print(dir," calculations started.")
    for n in range(1,n_repeats+1):
        #f = open(name+'/'+dir+f'/smd{n}_pullf.xvg','r')
        tf= loadtxt(name+'/'+dir+f'/smd{n}_pullf.xvg',skiprows=17)
        td= loadtxt(name+'/'+dir+f'/smd{n}_pullx.xvg',skiprows=17)




        if calculate_HB:
            ht = open(name + '/' + dir + f'/smd_hb_time{n}.dat', 'w')
            ht.write('# time [fs] vs number of hydrogen bonds between domains\n')
            # Hydrogen bonds number vs time
            #psf = [file for file in os.listdir(name) if file[-3:] == 'psf']
            print('Hydrogen bonds calculating - that may take some time...')
            u = Universe(name+'/'+pdb ,name+'/'+dir+f'/smd{n}.xtc')
            hbonds = HBA(universe=u,hydrogens_sel="protein and name H*",donors_sel=sel_const,acceptors_sel=sel_pull)
            hbonds.run()
            out = hbonds.count_by_time()
            hbonds_rev = HBA(universe=u,hydrogens_sel="protein and name H*",donors_sel=sel_pull,acceptors_sel=sel_const)
            hbonds_rev.run()
            out2 = hbonds_rev.count_by_time()
            t=0
            for i in (out+out2):
                ht.write(str(t)+' '+str(i)+'\n')
                t+= 0.05 #[ns]
            ht.close()



    #Writing down used pulling vectors
    theta = double(dir.split('_')[2])
    phi = double(dir.split('_')[4])
    x = cos(deg2rad(phi))*sin(deg2rad(theta))
    y = sin(deg2rad(phi))*sin(deg2rad(theta))
    z = cos(deg2rad(theta))
    bunch.append([x,y,z])
savetxt(name+'/bunch_of_vectors.dat',bunch)




# Plots - this script is adjusted for 9 pulling directions. If you want more, please increase the number of rows (9) and number of columns (3) accordingly
# Multiple repeats results are plotted as average value + outline based on calculated standard deviation value.

fig, ax = plt.subplots(3,3,sharex='col',figsize=(7,4))
#figv, axsv = plt.subplots(3, 1, sharex='col', figsize=(3, 4), subplot_kw={'projection': '3d'})
#ax[0,0] = fig.add_subplot( projection='3d')
i=0
sorted=['SMD_theta_0_phi_0','SMD_theta_45_phi_0','SMD_theta_45_phi_90','SMD_theta_45_phi_180','SMD_theta_45_phi_270','SMD_theta_90_phi_0','SMD_theta_90_phi_90','SMD_theta_90_phi_180','SMD_theta_90_phi_270'] #you might want to change order of your directions to plot

# looking for a maximum force value and HB numbers in all files to rescale figures
run = range(0,n_repeats)

max_force = max( [max([max(loadtxt(name+'/'+dir+f'/smd{n}_pullf.xvg',skiprows=17)[:,1]) for n in run]) for dir in sorted[:]] )
max_HB = max( [max([max(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')[:,1]) for n in run]) for dir in sorted[:]] )
print(max_force, max_HB)



scale = [250,640,389] #where your data stop (if you have different length of each repeat)
for dir, s in zip(sorted[:],scale):
    g = s*5
    gg = s
    print(dir)

    theta = double(dir.split('_')[2]) #these values are to draw bunch of vectors
    phi = double(dir.split('_')[4])



    for n in range(1,n_repeats+1):

        if n == 1:
            ft = reshape(loadtxt(name+'/'+dir+f'/smd{n}_pullf.xvg',skiprows=17)[:g,0],(len(loadtxt(name+'/'+dir+f'/smd{n}_pullf.xvg',skiprows=17)[:g]),1))
            ht = reshape(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')[:gg,0],(len(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')[:gg]),1))
            #fd = reshape(loadtxt(name+'/'+dir+f'/smd{n}_pullx.xvg',skiprows=17)[:,0],(len(loadtxt(name+'/'+dir+f'/smd{n}_pullf.xvg',skiprows=17)[:]),1))
            print(shape(ft),shape(ht))


        ft = concatenate((ft,reshape(loadtxt(name+'/'+dir+f'/smd{n}_pullf.xvg',skiprows=17)[:g,1],(len(loadtxt(name+'/'+dir+f'/smd{n}_pullf.xvg',skiprows=17)[:g]),1))),axis=1)
        ht = concatenate((ht, reshape(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')[:gg,1],(len(loadtxt(name+'/'+dir+f'/smd_hb_time{n}.dat')[:gg]),1))),axis=1)
        #fd = concatenate((fd, reshape(loadtxt(name+'/'+dir+f'/smd{n}_pullx.xvg',skiprows=17)[:,1],(len(loadtxt(name+'/'+dir+f'/smd{n}_pullf.xvg',skiprows=17)[:]),1))),axis=1)
        print(shape(ft),shape(ht))


    # filename = name+'/bunch_of_vectors.dat'  # Change to your filename
    # axsv[i].remove()
    # axsv[i] = figv.add_subplot(3,1,1+i,projection='3d')
    # plot_vectors(filename,theta,phi,axsv[i])




    print(shape(ft))
    for j in range(shape(ft)[0]):
        for k in [1,2,3]:
            if ft[j,k] < 0 or ft[j,k] > 600:
                ft[j,k] = ft[j-1,k]




    #Force vs time
    ax[i, 0].fill_between(ft[:g:red, 0] / 1000, mean(ft[:g:red,1:],axis=1) + std(ft[:g:red,1:],axis=1), mean(ft[:g:red,1:],axis=1) - std(ft[:g:red,1:],axis=1), linewidth=1,alpha=0.5)
    ax[i,0].plot(ft[:g:red,0]/1000,   mean(ft[:g:red,1:],axis=1)    ,linewidth=2)
    ax[i,0].set_xlim(0,10)
    ax[i, 0].set_ylim(0, 500)#max_force)
    #ax[i,0].set_xlabel('Time [ns]').set_fontsize(16)
    #ax[i,0].set_xticks(14)
    if i == 1: ax[i,0].set_ylabel('Force [pN]')
    #ax[i,0].yticks(size=14)


    #Force vs distance
    #ax[i, 1].fill_between(fd[:, 0], mean(fd[:, 1:], axis=1) + std(fd[:, 1:], axis=1),
    #                      mean(fd[:, 1:], axis=1) - std(fd[:, 1:], axis=1), linewidth=1, alpha=0.5)
    #ax[i,1].scatter(fd[::red,0],mean(fd[::red,1:],axis=1),s=2,c=ft[::red,0])
    #ax[i,1].set_xlim(30,55)
    #ax[i,1].set_ylim(0, max_force)
    #ax[i,1].set_xlabel(r'Distance [$\AA$]').set_fontsize(16)

    #if i == 4: ax[i,1].set_ylabel('Force [pN]')
    #

    # #hydrogen bonds number vs time

    ax[i, 2].fill_between(ht[:gg, 0]*0.2, mean(ht[:gg, 1:], axis=1) + std(ht[:gg, 1:], axis=1),
                          mean(ht[:gg, 1:], axis=1) - std(ht[:gg, 1:], axis=1), linewidth=1, alpha=0.5)
    ax[i,2].plot(ht[:gg,0]*0.2,mean(ht[:gg,1:],axis=1),linewidth=2)
    ax[i,2].set_xlim(0,10)
    ax[i,2].set_ylim(0,40)

    if i == 1: ax[i,2].set_ylabel('HB number')
    i+= 1


ax[2,0].set_xlabel('Time [ns]')
ax[2,1].set_xlabel(r'Distance [$\AA$]')
ax[2,2].set_xlabel('Time [ns]')

fig.tight_layout()
plt.savefig(name+'/all.png',dpi=600)


plt.show()
