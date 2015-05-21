import numpy as np
from astropy.io import ascii
from astropy import constants as const 
from astropy import units as u


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
print "file read"


# Now rescale all the columns to physical units
f40['T'] *= scale['tstar'].data[0]  # Myr
f40['MASS1'] *= scale['mstar'].data[0] # msun
f40['MASS2'] *= scale['mstar'].data[0] # msun
f40['RADIUS1'] *= scale['rstar'].data[0]*const.pc.cgs/const.R_sun.cgs # rsun
f40['RADIUS2'] *= scale['rstar'].data[0]*const.pc.cgs/const.R_sun.cgs # rsun
f40['A'] *= scale['rstar'].data[0]*const.pc.cgs/u.cm.cgs # cm
f40['P'] *= scale['tstar'].data[0] #Myr

mylist = ['T','NAME1','NAME2','KS1','KS2','MASS1','MASS2',
            'RADIUS1','RADIUS2','NP','E','A','P']

ascii.write(f40[mylist],'f40_physical.dat')
print "new file written"
