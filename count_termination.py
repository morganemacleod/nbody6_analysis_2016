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


#SCALING
TSTAR = 0.229

# CONTROL OPTIONS
plot_tres_hist = False



namelist57 = ('time','ID1','KS1','mass1','ID2','KS2','mass2','radius2','rt','rp','ecc')

print "starting to read 57 ..."
f57 = ascii.read("fort.57",
                 names=namelist57,
                 guess=False,
                 format='no_header',
                 delimiter=' ')
print "f57 read..."
#print f57
print "TDE rate [yr^-1] = ", len(f57['mass2'])/f57['time'][-1]/TSTAR*1.e-6  


namelist40 = ('T','NAME1','NAME2','KS1','KS2','MASS1','MASS2',
            'RADIUS1','RADIUS2','NP','E','A','P','EB','E3',
            'BODYI','CMX','CMY','CMZ','VCMX','VCMY','VCMZ'
            )

print "starting to read fort.40... "
f40 = ascii.read("fort.40",
                 names=namelist40,
                 guess=False,
                 format='no_header',
                 delimiter=' ')
print "file read"

#print f40


namelistinsp = ('time','ID1','KS1','mass1','ID2','KS2','mass2','radius2','rp','a','ecc')
insp = ascii.read("inspiral.dat",
                  names=namelistinsp,
                  guess=False,
                  format='no_header',
                  delimiter=' ')
print "inspirals read..."


# count the TDEs of MB
count_mb_TDEs = 0
for i,time in enumerate(f57['time']):
    # evaluate whether disrupted star is most-bound
    tfloor = np.floor(time)
    mb = f40[np.round(f40['T'])==tfloor]['NAME2']
    if (len(mb) == 0) : mb = [-99]
#    print len(mb)
#    print f57['ID2'][i], mb[0]
    if(f57['ID2'][i] == mb[0]) : 
        count_mb_TDEs = count_mb_TDEs + 1

count_TDEs = len(f57['time']) 
print "Out of %i TDEs, %i were the most bound, or a fraction: %f" % (count_TDEs,count_mb_TDEs,float(count_mb_TDEs)/float(count_TDEs))

print "There were %i compact object inspirals" % len(insp['time'])


namelist50 = ('nwrite','time','gamma',
            'prim_name','prim_mass','prim_ks','prim_radius',
            'comp_name','comp_mass','comp_ks','comp_radius',
            'A0','ECC','period',
            'pert_name','pert_mass','pert_ks','pert_radius',
            'peri_min','ECC1','VINF')

print "starting to read 50... "
f50 = ascii.read("fort.50",
                 names=namelist50,
                 guess=False,
                 format='no_header',
                 delimiter=' ')
print "fort.50 read..."
#print f50


# count the exchanges
count_ex=0
count_ex_tot=0 
t_residence = []
t_prev = 0
for i,time in enumerate(f50['time']) : 
    if(f50['prim_name'][i] == f50['prim_name'][i-1] and
       f50['comp_name'][i] != f50['comp_name'][i-1]) :
        dt = f50['time'][i] - t_prev
        t_residence = np.append(dt,t_residence)
        t_prev = f50['time'][i]
        count_ex_tot = count_ex_tot +1
        if(f50['comp_name'][i] == f50['pert_name'][i-1]) :
            count_ex = count_ex + 1
        
                        
print "This algorithm caught %i out of %i exchanges, or a fraction: %f" % (count_ex,count_ex_tot,1.0*count_ex/count_ex_tot) 
   

print t_residence

ascii.write([t_residence],"t_residence.dat",format='no_header')

if plot_tres_hist :
    plt.clf()
    plt.hist(np.log10(t_residence*TSTAR*1.e6),
             histtype='step',
             color='black',
             normed=True,
             bins=30)
    plt.xlabel(r'log companion residence time $[yr]$')
    plt.ylabel('PDF')
    plt.tight_layout(pad=1)
    plt.savefig('companion_t_residence_hist.eps',format='eps')    
