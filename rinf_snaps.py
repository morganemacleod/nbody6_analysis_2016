import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import ascii
import glob 

plot_save = True
plot_E_L = False
save_a_time = True
count_binarity = False


# set some global options
plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 16


#
# PARAMETERS
#
bhmassf = 0.002


#
# READ IN THE COMPANION FILE... 
#
times = ascii.read("snaps/times.dat",
                   format='no_header')
print times

namelist = ('T','NAME1','NAME2','KS1','KS2','MASS1','MASS2',
            'RADIUS1','RADIUS2','NP','E','A','P','EB','E3',
            'BODYI','CMX','CMY','CMZ','VCMX','VCMY','VCMZ'
            )

print "starting to read f40 file... "
f40 = ascii.read("fort.40",
                 names=namelist,
                 guess=False,
                 format='no_header',
                 delimiter=' ')
print "file read"






# IF READING THE FULL SNAPSHOTS AND SAVING THE INF RADII ONES
filenames = glob.glob('snaps/snap_0****_rinf.dat')
filenames = sorted(filenames)
print filenames
print ""
print "Reading ",len(filenames), "  files..."
print ""

#
# READ FILES
#
a_time = np.zeros((len(filenames),12))
comp_sma = np.zeros(len(filenames))
for j,filename in enumerate(filenames):
    print "reading ", filename
    snap = ascii.read(filename)
    Nparticles = len(snap)
    print "Particles in snapshot = ", Nparticles
    #print snap

    # get bh mass
    bhmass = np.interp(times['col1'][j],f40['T'],f40['MASS1'])
    print "bh mass = ",bhmass

    # get companion SMA 
    comp_sma[j] = np.interp(times['col1'][j],f40['T'],f40['A'])
    

    # Make Some Plots... 
    
    angular_momentum = np.zeros(Nparticles)
    ener = np.zeros(Nparticles)
    semi = np.zeros(Nparticles)
    for i in range(Nparticles) :

        pos = np.array([ snap['x'][i], snap['y'][i], snap['z'][i] ])
        vel = np.array([ snap['vx'][i],snap['vy'][i],snap['vz'][i] ])
        # particle angular momentum
        angular_momentum[i] =  np.sqrt(np.sum(np.cross(pos,vel)**2.))
        # particle energy
        ener[i] = (0.5 * np.sum(vel**2.) - bhmass / np.sqrt(np.sum(pos**2.)) )
        semi[i] = - bhmass /(2.0*ener[i])

        
    if plot_E_L :
        plt.plot(np.log10(angular_momentum),np.log10(semi),'.')
        plt.show()
        
    # Write the a-time array
    # offset to only include 2nd most bound...
    #if(min(snap['x']**2.) == 0.0) : 
    #    offset = 1
    #else :
    #    offset = 0
    offset = 1


    a_time[j][0] = times['col1'][j]
    a_time[j][1] = max( np.sqrt(snap['x']**2. + snap['y']**2. + snap['z']**2.))
    a_time[j][2:12] = sorted(semi[semi>0])[0+offset:10+offset]
    #print a_time[j][1],sorted(semi[semi>0])[0],sorted(semi[semi>0])[1] 
    #print sorted(semi[semi>0])
    #print "offset = ",offset, a_time[j][1],sorted(semi[semi>0])[0],sorted(semi[semi>0])[1] 

# now save / plot a vs time
if save_a_time:
    namesaout = ('time','rmbh','a1','a2','a3','a4','a5','a6','a7','a8','a9','a10')
    ascii.write(a_time,"semi_time_inner_10.dat",
                names=namesaout)


if count_binarity :
    sma_ratio = np.zeros(len(filenames))
    count_bin = 0
    for j in range(len(filenames)) :
        sma_ratio[j] = comp_sma[j]/a_time[j][2]
        if (sma_ratio[j] < 1./3. ) : count_bin = count_bin +1

    print ""
    print "BH is in a binary ",float(count_bin)/float(len(filenames))*100., " %  of the time"
    print ""


    print sma_ratio
    plt.hist(np.log10(sma_ratio),
             histtype='step',
             cumulative='True',
             normed=1,
             bins=30,
             color='black')
    plt.xlim(-3,0)
    plt.ylim(0,1)
    plt.xlabel(r"$\log \left(a_0/a_1 \right)$")
    plt.ylabel("CDF")
    plt.tight_layout(pad=1)
    #plt.show()
    plt.savefig("sma_hiarch_hist.eps")