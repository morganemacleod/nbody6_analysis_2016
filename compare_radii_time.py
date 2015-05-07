import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii


# set some global options
plt.rcParams['figure.figsize'] = (6,5)
plt.rcParams['legend.frameon'] = False
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['legend.labelspacing'] = 0.1
plt.rcParams['legend.handletextpad'] = 0.1
plt.rcParams['font.family'] = 'stixgeneral'
plt.rcParams['font.size'] = 16


folders = ["../king_w7rt_A/","../king_w7rt_B/","../king_w7rt_C/",
           "../king_w7_A/","../king_w7_B/","../king_w7_C/",
           "../king_w7lk_B/","../king_w7lk_C/","../king_w7lk_D/"]

chm_file = "core_hm_radii.dat"
time_file = "mttime.dat"

colors = ["DarkGreen","DarkGreen","DarkGreen",
          "Black","Black","Black",
          "DodgerBlue","DodgerBlue","DodgerBlue"]

for i,directory in enumerate(folders):
    
    print "reading from dir:  ",directory

    # READ MTTIME
    namelistmttime = ['mttime']
    mttime = ascii.read(directory+time_file,
                    names=namelistmttime,
                    guess=False,
                    format='no_header',
                    delimiter=' ')

    # READ CORE/HM RADII FILE
    namelistcorehm = ('core','hm')
    ch = ascii.read(directory+chm_file,
                    names=namelistcorehm,
                    guess=False,
                    format='no_header',
                    delimiter=' ')

    print "... files have length:", len(mttime), len(ch)
    
    # READ THE SCALES
    scales = ascii.read(directory+'scales_astropy.dat')
    tstar = scales['tstar'].data[0]
    rstar = scales['rstar'].data[0]


    plt.plot(mttime['mttime']*tstar/1.e3,ch['core']*rstar,'-',
             rasterized=True,
             linewidth=1.,
             color=colors[i],
             alpha=0.5)
    
    plt.plot(mttime['mttime']*tstar/1.e3,ch['hm']*rstar,'-',
             rasterized=True,
             linewidth=1.,
             color=colors[i],
             alpha=0.5)


plt.ylabel('Radius [pc]')
plt.xlabel('Time [Gyr]')
plt.tight_layout(pad=1)
plt.savefig('core_hm_comp.pdf')
