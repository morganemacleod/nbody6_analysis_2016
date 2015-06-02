import numpy as np
from astropy.io import ascii


scale = ascii.read('scales_astropy.dat')
print "scales read ... "

tr = ascii.read('t_residence.dat',format='no_header',names=['tr'])
print "tr read ... "
print tr

tr['tr'] *= scale['tstar'].data[0]

ascii.write(tr,'t_res_physical.dat',format='no_header')
print "physical scale file written"
