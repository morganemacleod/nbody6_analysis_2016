import numpy as np
from astropy.io import ascii
from astropy import constants as const 
from astropy import units as u
from sys import exit

scale = ascii.read('scales_astropy.dat')
print "scales read"

namelist = ('T','NAME1','NAME2','KS1','KS2','MASS1','MASS2',
            'RADIUS1','RADIUS2','NP','E','A','P','EB','E3',
            'BODYI','CMX','CMY','CMZ','VCMX','VCMY','VCMZ'
            )

print "starting to read file... "
f40 = ascii.read("fort.40",
                 names=namelist,
                 guess=False,
                 format='no_header',
                 delimiter=' ')
print "... fort.40 read"


# read in the evolve file and prepare to add that data
print "starting to read sse database..."
sse = ascii.read("/trove/mmacleod/nbody_runs/nbody_analysis/evolve_mo.dat",
                 delimiter=' ',
                 guess=False)
print "... sse read"

# initial mass in msun, age in Myr
def LRT_at_age(Mi,age,decimal=1):
    mask = (sse['Mi']==np.round(Mi,decimals=decimal))
    Lum = np.interp(age,sse[mask]['Tev'],sse[mask]['logL'])
    Rad = np.interp(age,sse[mask]['Tev'],sse[mask]['logR'])
    Teff = np.interp(age,sse[mask]['Tev'],sse[mask]['logT'])

    return Lum,Rad,Teff


# read in the zeroth snapshot
namelistsnap = ('mass','x','y','z','vx','vy','vz','ID','KS','radius')
print "about to read snap"
snap = ascii.read("snaps/snap_00001.dat",
                  guess=False,
                  format='no_header',
                  names=namelistsnap,
                  delimiter=' ')
print "...snap read"

# initial mass in NB units
def mo(ID):
    mask = (snap['ID']==ID)
    return snap[mask]['mass'].data[0]



# Now rescale all the columns to physical units
f40['T'] *= scale['tstar'].data[0]  # Myr
f40['MASS1'] *= scale['mstar'].data[0] # msun
f40['MASS2'] *= scale['mstar'].data[0] # msun
f40['RADIUS1'] *= scale['rstar'].data[0]*const.pc.cgs/const.R_sun.cgs # rsun
f40['RADIUS2'] *= scale['rstar'].data[0]*const.pc.cgs/const.R_sun.cgs # rsun
f40['A'] *= scale['rstar'].data[0]*const.pc.cgs/u.cm.cgs # cm
f40['P'] *= scale['tstar'].data[0] #Myr


# fill in some new columns
f40['MASS20'] = 0.0
f40['RAD2'] = 0.0
f40['LUM2'] = 0.0
f40['TEFF2'] = 0.0
for i in range(len(f40)):
    f40['MASS20'][i] = mo(f40['NAME2'][i])*scale['mstar'].data[0]
    L,R,T = LRT_at_age(f40['MASS20'][i],f40['T'][i])
    f40['LUM2'][i] = 10**L
    f40['RAD2'][i] = 10**R
    f40['TEFF2'][i] = 10**T



mylist = ['T','NAME1','NAME2','KS1','KS2','MASS1','MASS2',
          'RADIUS1','RADIUS2','NP','E','A','P','MASS20',
          'LUM2','RAD2','TEFF2']

ascii.write(f40[mylist],'f40_physical.dat')
print "new file written"
