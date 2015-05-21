# Program to Read snapshots grabbed from snapdata.dat
# Prerequisites:
# mkdir snaps
# python split_snaps.py   -- creates the individual series of snaps
# grep "found at" output.dat | cut -c 16- > bh.dat   -- create bh.dat

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import ascii
import glob 
#from scipy.interpolate import interp1d

plot_cm_pos = False
plot_bh_rel_cm = False

calc_profiles = True

plot_menc_radius = False
plot_menc_radius_bin = False

plot_dens_radius = False
plot_dens_radius_bin = False
plot_dens_radius_com_bh_comp = False

plot_ndens_radius = False
plot_ndens_radius_bin = False

plot_mavg_radius = False

plot_bh_ener = False
plot_save = False

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
#
npart0 = 100001
bhmassf = 0.001
# SCALINGS 
scale = ascii.read("scales_astropy.dat")
RSTAR = scale['rstar'].data[0]      # RADIUS OF 1 NB UNIT IN PC
MSTAR = scale['mstar'].data[0]  # TOTAL MASS, initial
VSTAR = scale['vstar'].data[0]  # velocity in KM/S
MAVG = scale['mavg'].data[0]  # Avg mass in msun
SU = scale['su'].data[0]   # NB length to solar radii ? 
TSTAR = scale['tstar'].data[0] # TIMESCALE IN MYR


filenames = glob.glob('snaps/snap_0[0-9][0-9][0-9][0-9].dat')
filenames = sorted(filenames)
print filenames
print ""
print "Reading ",len(filenames), "  files..."
print ""
# Setup some array-binning 
binwidth = 1
nbins = int(len(filenames)/binwidth)
print ""
print "binning ", len(filenames)," files into ",nbins, " bins of width ",binwidth 
print ""


#
# READ FILES
#
namelist = ('mass','x','y','z','vx','vy','vz','ID','KS','radius')


# Create some arrays
cm_data = np.zeros((len(filenames),4))
vcm_data = np.zeros((len(filenames),4))
bh_pos =  np.zeros((len(filenames),4))
bh_offset = np.zeros((len(filenames),4))
bh_vrel = np.zeros((len(filenames),4))
r_inf = np.zeros((len(filenames),2))

# define some data
if calc_profiles :
    mpoints = 120
    menc = np.zeros((len(filenames),mpoints))
    menc_bh = np.zeros((len(filenames),mpoints))
    ndens = np.zeros((len(filenames),mpoints))
    ndens_bh = np.zeros((len(filenames),mpoints))
    dens = np.zeros((len(filenames),mpoints))
    dens_bh = np.zeros((len(filenames),mpoints))
    

    
# loop through reading files
for i,filename in enumerate(filenames):
    #print "FILE {}, {}".format(i,filename)
    #print "starting to read file... "
    snap = ascii.read(filename,
                      names=namelist,
                      guess=False,
                      format='no_header',
                      delimiter=' ')
    #print "file read"
    #print snap

    # Exclude the BH
    #snap_nobh = snap[:][snap['ID']>1]
    #snap_bh = snap[:][snap['ID']==1]
    snap_bh = snap[:][snap['mass']>0.99*bhmassf]
    snap_bh = snap_bh[:][snap_bh['KS']==14]
    print snap_bh
    #snap_bh = snap_bh[:][snap_bh['ID']<=npart0]
    snap_nobh = snap[:][snap['mass']<bhmassf]
    #print snap[:][snap['ID']>=npart0]
        
    # center of mass of the stars
    cm = [np.sum(snap_nobh['x']*snap_nobh['mass']),
          np.sum(snap_nobh['y']*snap_nobh['mass']),
          np.sum(snap_nobh['z']*snap_nobh['mass'])]/np.sum(snap_nobh['mass'])
    # velocity of center of mass of the stars
    vcm = [np.sum(snap_nobh['vx']*snap_nobh['mass']),
          np.sum(snap_nobh['vy']*snap_nobh['mass']),
          np.sum(snap_nobh['vz']*snap_nobh['mass'])]/np.sum(snap_nobh['mass'])
    
    
    cm_data[i,:] = [i,cm[0],cm[1],cm[2]]
    vcm_data[i,:] = [i,vcm[0],vcm[1],vcm[2]]
    bh_pos[i,:] = [i,snap_bh['x'],snap_bh['y'],snap_bh['z']]
    bh_offset[i,:] = [i,snap_bh['x']-cm[0],snap_bh['y']-cm[1],snap_bh['z']-cm[2] ]
    #bh_vrel[i,:] = [i,snap_bh['x']]
    print "cm = ", cm_data[i]
    #print "bh = ", bh_pos[i]

    
    if calc_profiles:

    
        #radii = np.linspace(1.0/mpoints,2.0,num=mpoints)
        radii = np.logspace(-3.,1.5,num=mpoints)
        #print radii
        for j,rad in enumerate(radii) :
            # set the inner radius
            rad_in = radii[j-1]
            # handle the first-index case
            if(j==0) : rad_in = 0.0
            # sum the mass of objects within the radius
            mask_out = (np.sqrt(
                    (snap_nobh['x']-cm[0])**2. + 
                    (snap_nobh['y']-cm[1])**2. + 
                    (snap_nobh['z']-cm[2])**2.)<rad) 
            mask_in = (np.sqrt((snap_nobh['x']-cm[0])**2. + 
                              (snap_nobh['y']-cm[1])**2. + 
                              (snap_nobh['z']-cm[2])**2.)>=rad_in) 
            # relative to the BH
            mask_out_bh = (np.sqrt(
                    (snap_nobh['x']-snap_bh['x'])**2. + 
                    (snap_nobh['y']-snap_bh['y'])**2. + 
                    (snap_nobh['z']-snap_bh['z'])**2.)<rad) 
            mask_in_bh = (np.sqrt((snap_nobh['x']-snap_bh['x'])**2. + 
                                  (snap_nobh['y']-snap_bh['y'])**2. + 
                                  (snap_nobh['z']-snap_bh['z'])**2.)>=rad_in) 
            # select everything within the radius (solar masses)
            selection_out = snap_nobh['mass'][mask_out] * MSTAR
            menc[i,j] = np.sum( selection_out )
            selection_out_bh = snap_nobh['mass'][mask_out_bh] * MSTAR
            menc_bh[i,j] = np.sum( selection_out_bh )
            # volume element (pc^3)
            dVol = 4./3. * np.pi * (rad**3.0 - rad_in**3.0 ) * RSTAR**3
            # select between radii (solar masses) 
            selection = snap_nobh["mass"][mask_out & mask_in] * MSTAR
            selection_bh = snap_nobh["mass"][mask_out_bh & mask_in_bh] * MSTAR 
            # com case 
            Mstars = np.sum(selection)
            Nstars = len(selection)
            # bh case
            Mstars_bh = np.sum(selection_bh)
            Nstars_bh = len(selection_bh)
            # allocate ndens and dens
            dens[i,j] = Mstars/dVol
            ndens[i,j] = Nstars/dVol
            # bh case
            dens_bh[i,j] = Mstars_bh/dVol
            ndens_bh[i,j] = Nstars_bh/dVol
            
        #
        # SPHERE OF INF DIAGNOSTICS
        #
        # now that have calculated the profile, compute the menc=mbh sphere of influence
        radius_of_influence  = np.interp(snap_bh['mass'][0]*MSTAR,menc_bh[i],radii)
        r_inf[i,:] = [i,radius_of_influence]
        print "r_inf = ",r_inf[i]
        mask_rinf = (np.sqrt(
                (snap_nobh['x']-snap_bh['x'])**2. +
                (snap_nobh['y']-snap_bh['y'])**2. +
                (snap_nobh['z']-snap_bh['z'])**2.)<r_inf[i,1])
        selection_rinf = snap_nobh[mask_rinf]
        print "N* (<rinf) = ", len(selection_rinf)
        print "M* (<rinf) = ", np.sum(selection_rinf["mass"]*MSTAR) 
        #print "next most massive stars =", sorted(selection_rinf['mass']*MSTAR)[-1:-3]
        
        # rank the energies
        energies_rinf = (0.5 * (selection_rinf['vx']**2. + 
                                selection_rinf['vy']**2. + 
                                selection_rinf['vz']**2.) -
                         snap_bh['mass'][0] / np.sqrt((selection_rinf['x']-snap_bh['x'][0])**2.+
                                                      (selection_rinf['y']-snap_bh['y'][0])**2.+
                                                      (selection_rinf['z']-snap_bh['z'][0])**2.)
                         )
        print "radii (except companion) = ", sorted( np.sqrt((selection_rinf['x']-snap_bh['x'][0])**2.+
                                                             (selection_rinf['y']-snap_bh['y'][0])**2.+
                                                             (selection_rinf['z']-snap_bh['z'][0])**2.) )[1:5]
        #print "dx min = ", sorted( abs(selection_rinf['x']-snap_bh['x'][0]) )[0]
        #print "dy min = ", sorted( abs(selection_rinf['y']-snap_bh['y'][0]) )[0]
        #print "dz min = ", sorted( abs(selection_rinf['z']-snap_bh['z'][0]) )[0]



        print "BE (except companion) = ", sorted(energies_rinf)[1:6]
        mb_ind = np.argmin(energies_rinf)
        print "most bound = ",selection_rinf["ID"][mb_ind],selection_rinf["mass"][mb_ind] 
        #print "X =", selection_rinf['x'][mb_ind], snap_bh['x'][0]
        #print "Y =", selection_rinf['y'][mb_ind], snap_bh['y'][0]
        #print "Z =", selection_rinf['z'][mb_ind], snap_bh['z'][0]
        #print ""
        print ""
        # for saving the snapshot, modify the positions & velocities
        selection_rinf['x'] = selection_rinf['x'] - snap_bh['x'][0]
        selection_rinf['y'] = selection_rinf['y'] - snap_bh['y'][0]
        selection_rinf['z'] = selection_rinf['z'] - snap_bh['z'][0]
        selection_rinf['vx'] = selection_rinf['vx'] - snap_bh['vx'][0]
        selection_rinf['vy'] = selection_rinf['vy'] - snap_bh['vy'][0]
        selection_rinf['vz'] = selection_rinf['vz'] - snap_bh['vz'][0]
        rinf_filename = filename[0:-4]+"_rinf.dat"
        ascii.write(selection_rinf,rinf_filename)
        print "saved ",rinf_filename

        #print menc[i]
        
        



# binned arrays
menc_bin = np.zeros((nbins,mpoints))
menc_bin_bh = np.zeros((nbins,mpoints))
dens_bin = np.zeros((nbins,mpoints))
dens_bin_bh = np.zeros((nbins,mpoints))
ndens_bin = np.zeros((nbins,mpoints))
ndens_bin_bh = np.zeros((nbins,mpoints))

# bin slices
slices = np.linspace(0,len(filenames),nbins,endpoint=False).astype(np.int)
print slices

for j,rad in enumerate(radii) :
    menc_bin[:,j] = np.add.reduceat(menc[:,j],slices)/binwidth
    dens_bin[:,j] = np.add.reduceat(dens[:,j],slices)/binwidth
    ndens_bin[:,j] = np.add.reduceat(ndens[:,j],slices)/binwidth
    
    # bh case
    menc_bin_bh[:,j] = np.add.reduceat(menc_bh[:,j],slices)/binwidth    
    dens_bin_bh[:,j] = np.add.reduceat(dens_bh[:,j],slices)/binwidth
    ndens_bin_bh[:,j] = np.add.reduceat(ndens_bh[:,j],slices)/binwidth


# Read in bh-track file
# grep "found at" output.dat | cut -c 16- > bh.dat
#namelistbh = ('ID','Mass','X','Y','Z')
#bh = ascii.read("bh.dat",
#                names=namelistbh)
#print "BH file read"

print "reading times ..."
namelisttime = ('time')
times = ascii.read("snaps/times.dat",format='no_header') 
#print times


#print bh_pos


if plot_cm_pos:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(cm_data[:,1],cm_data[:,2],cm_data[:,3])
    #ax.plot(bh_pos[:,1],bh_pos[:,2],bh_pos[:,3])
    ax.plot(bh['X'],bh['Y'],bh['Z'])
    plt.show()


if plot_bh_rel_cm:
    
    
    #cmxf = interp1d(times['col1'], cm_data[:,1], kind='cubic')
    #cmyf = interp1d(times['col1'], cm_data[:,2], kind='cubic')
    #cmzf = interp1d(times['col1'], cm_data[:,3], kind='cubic')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(bh_offset[:,1],bh_offset[:,2],bh_offset[:,3],'.')
    #ax.plot(bh_pos[:,1]-cm_data[:,1],bh_pos[:,2]-cm_data[:,2],bh_pos[:,3]-cm_data[:,3])
    
    ax.set_xlim(-0.15,0.15)
    ax.set_ylim(-0.15,0.15)
    ax.set_zlim(-0.15,0.15)
    plt.show()


    plt.clf()
    plt.hist( np.log10(np.sqrt(bh_offset[:,1]**2 + bh_offset[:,2]**2 + bh_offset[:,3]**2)),
              normed=True,
              histtype='step')
    plt.show()

# SET COLORMAP FOR ALL PROFILE PLOTS
profile_cm = plt.cm.jet



if plot_menc_radius :
    plt.clf()
    for i in range(len(filenames)):
        plt.loglog(radii*RSTAR,menc_bh[i,:],
                   color=profile_cm(float(i)/len(filenames)) 
                   )
    plt.xlabel(r"$ r {\rm [pc]}$")
    plt.ylabel(r"$ M_{\rm enc} [M_\odot]$")
    plt.tight_layout(pad=1)
    #plt.show()
    plt.savefig("menc_radius.eps")

if plot_menc_radius_bin :
    plt.clf()
    for i in range(nbins):
        plt.loglog(radii*RSTAR,menc_bin_bh[i,:],
                   color=profile_cm(float(i)/len(filenames)) 
                   )
    plt.xlabel(r"$ r {\rm [pc]}$")
    plt.ylabel(r"$ M_{\rm enc} [M_\odot]$")
    plt.tight_layout(pad=1)
    #plt.show()
    plt.savefig("menc_radius_bin.eps")

if plot_dens_radius :
    plt.clf()
    for i in range(len(filenames)):
        plt.loglog(radii*RSTAR,dens_bh[i,:],
                   color=profile_cm(float(i)/len(filenames)) 
                   )

    plt.ylabel(r"$ \rho [M_\odot {\rm pc}^{-3}]$")
    plt.xlabel(r"$ r {\rm [pc]}$")    
    plt.tight_layout(pad=1)
    plt.show()

if plot_dens_radius_bin :
    plt.clf()
    for i in range(nbins):
        plt.loglog(radii*RSTAR,dens_bin_bh[i,:],
                   color=profile_cm(float(i)/len(filenames)) 
                   )
    plt.ylabel(r"$ \rho [M_\odot {\rm pc}^{-3}]$")
    plt.xlabel(r"$ r {\rm [pc]}$")    
    plt.tight_layout(pad=1)
    plt.savefig("dens_radius_bin.eps")

if plot_ndens_radius :
    plt.clf()
    for i in range(len(filenames)):
        plt.loglog(radii*RSTAR,ndens_bh[i,:],
                   color=profile_cm(float(i)/len(filenames)) 
                   )
    plt.ylabel(r"$ n [{\rm pc}^{-3}]$")
    plt.xlabel(r"$ r {\rm [pc]}$")
    plt.tight_layout(pad=1)
    #plt.show()
    plt.savefig("ndens_radius.eps")

if plot_ndens_radius_bin :
    plt.clf()
    for i in range(nbins):
        plt.loglog(radii*RSTAR,ndens_bin_bh[i,:],
                   color=profile_cm(float(i)/len(filenames)) 
                   )
    plt.xlabel(r"$ r {\rm [pc]}$")
    plt.ylabel(r"$ n [{\rm pc}^{-3}]$")
    plt.tight_layout(pad=1)
    #plt.show()
    plt.savefig("ndens_radius_bin.eps")

if plot_mavg_radius :
    plt.clf()
    for i in range(nbins):
        plt.loglog(radii*RSTAR,dens_bin_bh[i,:]/ndens_bin_bh[i,:],
                   color=plt.cm.rainbow(float(i)/len(filenames)) 
                   )
    plt.ylabel(r"$ \bar M [M_\odot]$")
    plt.xlabel(r"$ r {\rm [pc]}$")
    plt.tight_layout(pad=1)
    #plt.show()
    plt.savefig("mavg_radius_bin.eps")

if plot_dens_radius_com_bh_comp :
    plt.clf()
    for i in range(nbins):
        plt.loglog(radii*RSTAR,dens_bin[i,:],
                   color=profile_cm(float(i)/len(filenames)) 
                   )
        plt.loglog(radii*RSTAR,dens_bin_bh[i,:],
                   color='r'
                   )
    plt.ylabel(r"$ \rho [M_\odot {\rm pc}^{-3}]$")
    plt.xlabel(r"$ r {\rm [pc]}$")
    plt.tight_layout(pad=1)
    plt.show()

if plot_bh_ener :
    plt.clf()
    menc_rbh = np.zeros(len(filenames))
    for i in range(len(filenames)):
        enclosed_mass = interp1d(radii,menc[i,:])
        r_bh = np.sqrt(bh_offset[i,1]**2 + bh_offset[i,2]**2 + bh_offset[i,3]**2)
        menc_rbh[i] = enclosed_mass(r_bh)
        print menc_rbh[i]

    plt.hist(np.log10(menc_rbh[menc_rbh>0]/0.01),normed=True,histtype='step')
    plt.show()
