# Program to read HR-diagram diagnostics grabbed from fort.83
# to use:
# mkdir hr
# python split_hr.py

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import glob 
#from scipy.interpolate import interp1d

plot_mass_hist = True
plot_hr_series = True

filenames = glob.glob('hr/hr_single_00***.dat')

#
# Set up plots
#
if plot_mass_hist :
    mass_hist = plt.figure()
    mass_hist_ax = mass_hist.add_subplot(111)



# READ FILES
namelist = ('name','kstar','RI','mass','logL','logR','logTeff')
filenames = sorted(filenames)

for i,filename in enumerate(filenames):
    # read the file
    hr =  ascii.read(filename,
                     names=namelist,
                     guess=False,
                     format='no_header',
                     delimiter=' ')
    
    #print hr

    # PLOT A MASS HISTOGRAM WITH TIME
    if plot_mass_hist :
        cmap = plt.cm.RdYlBu_r
        mmin = 0.2
        mmax = 15.0
        mbins = np.linspace(np.log10(mmin),np.log10(mmax),15)
        mass_hist_ax.hist(np.log10(hr['mass']),
                          normed=False,
                          log=True,
                          histtype='step',
                          color=cmap( i/float(len(filenames)) ),
                          bins=mbins,
                          linewidth=2
                          )
        mass_hist_ax.set_xlabel("log mass [msun]")
        #print cmap( i/float(len(filenames)) )


    # PLOT A SERIES OF HR-DIAGRAMS
    if plot_hr_series :
        print "plotting figure ... %05i" % (i+1)
        hrfig = plt.figure(figsize=(5,7))
        hrax = hrfig.add_subplot(111)
        hrax.plot(hr['logTeff'],hr['logL'],'.')
        hrax.set_xlim(3.0,5.0)
        hrax.set_ylim(-4.0,5.0)
        hrax.set_xlabel("log Teff")
        hrax.set_ylabel("log L")
        hrfig.savefig("hr/hr_%05i.png" % (i+1),format='png')
        plt.close(hrfig)


# show/save the plots
if plot_mass_hist:
    mass_hist.savefig("mass_histogram_time.eps",format='eps')
    
