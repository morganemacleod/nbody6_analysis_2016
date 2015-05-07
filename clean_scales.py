import re
from astropy.io import ascii
import numpy as np

fin = 'scales.dat'
fout = 'scales_astropy.dat'
folders = ["../king_w7rt_A/","../king_w7rt_B/","../king_w7rt_C/",
           "../king_w7_A/","../king_w7_B/","../king_w7_C/",
           "../king_w7lk_B/","../king_w7lk_C/","../king_w7lk_D/"]

numbers = re.compile("-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-+]?\ *[0-9]+)?")

namelist = ['rstar','mstar','vstar','tstar','mavg','su']

for i,directory in enumerate(folders):
    
    print "reading from dir:  ",directory
    f = open(directory+fin,'r')
    char = f.readline()
    scalelist = np.array([numbers.findall(char)])
    print scalelist
    
    ascii.write(scalelist,directory+fout,names=namelist)
