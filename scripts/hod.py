
import caesar
import matplotlib as mpl
mpl.use('Agg')
import matplotlib
# Specify renderer
# matplotlib.use('Agg')
# Boiler-plate settings for producing pub-quality figures
# 1 point = 1/72 inch
matplotlib.rcParams.update({'figure.figsize'  : (8,5)    # inches
                            ,'font.size'       : 22      # points
                            ,'legend.fontsize' : 16      # points
                            ,'lines.linewidth' : 2       # points
                            ,'axes.linewidth'  : 2       # points
                            ,'text.usetex'     : True    # Use LaTeX to layout text
                            ,'font.family'     : "serif" # Use serifed fonts
                            ,'xtick.major.size'  : 6     # length, points
                            ,'xtick.major.width' : 2     # points
                            ,'xtick.minor.size'  : 3     # length, points
                            ,'xtick.minor.width' : 1     # points
                            ,'ytick.major.size'  : 6     # length, points
                            ,'ytick.major.width' : 2     # points
                            ,'ytick.minor.size'  : 3     # length, points
                            ,'ytick.minor.width' : 1     # points
                            ,'font.serif'      : ("Times"
                                                  ,"Palatino"
                                                  ,"Computer Modern Roman"
                                                  ,"New Century Schoolbook"
                                                  ,"Bookman")
                            ,'font.sans-serif' : ("Helvetica"
                                                  ,"Avant Garde"
                                                  ,"Computer Modern Sans serif")
                            ,'font.monospace'  : ("Courier"
                                                  ,"Computer Modern Typewriter")
                            ,'font.cursive'    : "Zapf Chancery"
})

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
import h5py
import sys
from astropy.cosmology import FlatLambdaCDM
from scipy import stats


inoutName = {
            'm50n512_036': '/mnt/ceph/users/daisyleung/simba/sim/m50n512/s50/Groups/m50n512_036.hdf5',
            'm100n1024_036': '/mnt/ceph/users/daisyleung/simba/sim/m100n1024/s50/Groups/m100n1024_036.hdf5',
            'm50n1024_036': '/mnt/ceph/users/daisyleung/simba/sim/m50n1024/s50/Groups/m50n1024_036.hdf5',
            'm25n1024_036': '/mnt/ceph/users/daisyleung/simba/sim/m25n1024/s50_new/Groups/m25n1024_036.hdf5',
            'm25n256_036': '/mnt/ceph/users/daisyleung/simba/sim/m25n256/s50/Groups/m25n256_036.hdf5'
            }

nbin = 20
allhalos = True    # False

for k, infile in inoutName.items():
    figname = k + '_hod_bin' + str(nbin) + '.pdf'

    sim = caesar.load(infile) # load caesar file
    mlim = 64*sim.simulation.critical_density.value*sim.simulation.boxsize.value**3*sim.simulation.omega_matter/sim.simulation.effective_resolution**3 # galaxy mass resolution limit in DM particle

    vol_mpc = sim.simulation.boxsize.to('Mpccm').d**3
    cosmo = FlatLambdaCDM(H0=100*sim.simulation.hubble_constant, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)

    myobjs = sim.galaxies
    cents = np.asarray([i for i in myobjs if i.central==1])
    Mhmin = np.min( np.asarray([i.halo.masses['total'] for i in cents]) ) # smallest halo to host a galaxy in each cental halo
    if not allhalos:
        # only include halos (central + satellites) that are above the min. central halo mass
        halos = np.asarray([i for i in sim.halos if i.masses['total']>Mhmin])
    else:
        figname = k + '_hod_bin' + str(nbin) + '_allhalos.pdf'
        halos = np.asarray([i for i in sim.halos])

    mh = np.asarray([i.masses['total'] for i in halos])

    pos = np.asarray([i.pos for i in halos])
    print('Mhalo,min=',np.log10(Mhmin),' Ncents=',len(cents),' Nhalos=',len(halos))

    ms = np.asarray([i.masses['stellar'] for i in myobjs])
    logms = np.log10(ms)
    sfr = np.asarray([i.sfr for i in myobjs])
    ssfr = np.log10(1.e9*sfr/ms+10**(-2.9+0.3*sim.simulation.redshift)) # with a floor to prevent NaN's

    halogals = np.asarray([i.galaxy_index_list for i in halos])
    halomags = []
    for i in range(len(halogals)):
        halomags.append(halogals[i])
    ngal = np.asarray([len(hm) for hm in halomags])


    # plot
    plt.close('all')
    plt.figure()
    import plotmedian as pm
    # I just estimate it based on the spread on individual quantities over 8 octants (defined inside plotmedian).  8 is purely for convenience, and seemed like a sufficiently large number to get a decent variance.  I’ve never compared to calculating it from the correlation function, which would probably be more accurate.
    pm.plotmedian(np.log10(mh),np.log10(ngal),ltype='-',lw=2,
                    ax=plt,
                    bins=nbin,
                    stat='mean',
                    pos=pos,
                    boxsize=sim.simulation.boxsize,
                    label=r'$z$ = %g'%(np.round(sim.simulation.redshift,1))
                  )

    plt.minorticks_on()
    # plt.xlim(11.0,)
    plt.ylabel(r'$<$log N$>$',fontsize=16)
    plt.xlabel(r'log $M_{\rm halo}$ [$M_{\odot}$]' ,fontsize=16)
    plt.legend(loc='upper left')

    plt.savefig(figname, bbox_inches='tight', format='pdf')
    plt.show()
