"""

Determine cosmic variance (as a func of survey area), from TPCF (integrating over relevant columns)


Python: http://www.mpia.de/homes/dgoulier/MLClasses/Astro-Tutorial%20-%20The%20Two-Point%20Correlation%20Function%20with%20R%20and%20Python.html#the_two-point_correlation_function

https://halotools.readthedocs.io/en/latest/quickstart_and_tutorials/tutorials/catalog_analysis/galcat_analysis/basic_examples/clustering_examples/galaxy_catalog_analysis_tutorial2.html#galaxy-catalog-analysis-tutorial2


read Coil13_LSS.pdf


"""

import numpy as np
import pandas as pd
import numpy as np
import matplotlib

# Specify renderer
matplotlib.use('Agg')

# Boiler-plate settings for producing pub-quality figures
# 1 point = 1/72 inch
from cycler import cycler
matplotlib.rcParams['axes.prop_cycle'] = cycler(color='bgrcmyk')

matplotlib.rcParams.update({'figure.figsize': (8, 5)    # inches
                            , 'font.size': 22      # points
                            , 'legend.fontsize': 16      # points
                            , 'lines.linewidth': 1.5       # points
                            , 'axes.linewidth': 1.5       # points
                            , 'text.usetex': True    # Use LaTeX to layout text
                            , 'font.family': "serif"  # Use serifed fonts
                            , 'xtick.major.size': 10     # length, points
                            , 'xtick.major.width': 1.5     # points
                            , 'xtick.minor.size': 6     # length, points
                            , 'xtick.minor.width': 1     # points
                            , 'ytick.major.size': 10     # length, points
                            , 'ytick.major.width': 1.5     # points
                            , 'ytick.minor.size': 6     # length, points
                            , "xtick.minor.visible": True
                            , "ytick.minor.visible": True
                            , 'font.weight': 'bold'
                            , 'ytick.minor.width': 1     # points
                            , 'font.serif': ("Times", "Palatino", "Computer Modern Roman", "New Century Schoolbook", "Bookman"), 'font.sans-serif': ("Helvetica", "Avant Garde", "Computer Modern Sans serif"), 'font.monospace': ("Courier", "Computer Modern Typewriter"), 'font.cursive': "Zapf Chancery"
                            })


import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

# read in data
import caesar
nn = 'm50n1024'
# nn = 'm100n1024'
# nn = 'm25n1024'

infile = '/mnt/ceph/users/daisyleung/simba/sim/' + nn + '/s50/Groups/' + nn + '_036.hdf5'
sim = caesar.load(infile)  # load caesar file
mlim = 64 * sim.simulation.critical_density.value * sim.simulation.boxsize.value**3 * \
    sim.simulation.omega_matter / \
    sim.simulation.effective_resolution**3  # galaxy mass resolution limit in DM particle

boxsize = sim.simulation.boxsize.to('Mpccm').d
vol_mpc = boxsize**3
cosmo = FlatLambdaCDM(H0=100 * sim.simulation.hubble_constant,
                      Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)

# # =========================================================
# import numpy as np
# from os.path import dirname, abspath, join as pjoin
# from Corrfunc.theory.DD import DD
# from Corrfunc.utils import convert_3d_counts_to_cf

# # Generate randoms on the box
# N = len(galpos)
# rand_N = 3*N
# rand_X = np.random.uniform(0, boxsize, rand_N)
# rand_Y = np.random.uniform(0, boxsize, rand_N)
# rand_Z = np.random.uniform(0, boxsize, rand_N)


# # Auto pair counts in DD
# nbins = 10
# bins = np.linspace(0.1, 10.0, nbins + 1) # note that +1 to nbins
# nthreads= 10

# autocorr=1
# DD_counts = DD(autocorr, 10, bins, galpos[:, 0]/1.e3, galpos[:, 1]/1.e3,
#                galpos[: 2]/1.e3,
#                periodic=False, verbose=True)


# # Cross pair counts in DR
# autocorr=0
# DR_counts = DD(autocorr, 10, bins, galpos[:, 0], galpos[:, 1]/1.e3,
#                galpos[: 2],
#                X2=rand_X, Y2=rand_Y, Z2=rand_Z,
#                periodic=False, verbose=True)

# # Auto pairs counts in RR
# autocorr=1
# RR_counts = DD(autocorr, 10, bins, rand_X, rand_Y, rand_Z,
#                periodic=False, verbose=True)

# # All the pair counts are done, get the correlation function
# cf = convert_3d_counts_to_cf(N, N, rand_N, rand_N,
#                             DD_counts, DR_counts,
#                             DR_counts, RR_counts)



# --------------- halotools ---------------
import halotools
from halotools.mock_observables import return_xyz_formatted_array
from halotools.mock_observables import tpcf, tpcf_jackknife
from halotools.mock_observables import wp, rp_pi_tpcf
import numpy as np
from halotools.mock_observables import tpcf_one_two_halo_decomp


# TPCF (r, m, z), where z fixed at z=6, we take m as the median of all gal
# sample centered on what stellar mass range
plt.figure()
mstar = np.array([g.masses['stellar'] for g in sim.galaxies])
plt.hist(mstar, bins=100)
plt.title('median: {:.1f}E8 Msun; delta: {:.1f}E10 Msun'.format(np.median(mstar)/1.e8, (mstar.max() - mstar.min())/1.e10))
plt.xscale('log')
plt.yscale('log')
plt.savefig('mstar_dist' + nn + '.pdf', bbox_inches='tight')
# import sys; sys.exit()

# ACF
galpos = np.array([g.pos.d for g in sim.galaxies])/1.e3      # cMpc/h
#  the real space radial bins in which pairs are counted
rsize = 15
rbins = np.logspace(-1, 1.25, rsize)                      # cMpc/h
rbin_centers = (rbins[1:] + rbins[:-1]) / 2.
# real space TPCF
xi_all = tpcf(galpos, rbins, period=boxsize,
              num_threads='max', estimator='Landy-Szalay')

plt.figure()
plt.plot(rbin_centers, xi_all,
         label='All galaxies', color='k', marker='o')
plt.loglog()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r'$r $  $\rm{[Mpc h^{-1}]}$', fontsize=25)
plt.ylabel(r'$\\xi(r)$', fontsize=25)
plt.legend(loc='best', fontsize=20)
plt.savefig('realspace_xi_' + nn + '.pdf', bbox_inches='tight')


# fit power law to correlation
from scipy.optimize import curve_fit

def TPCF_powerlaw(r, c, r0, gamma):
    """
    r0:
        correlation length

    gamma:
        power law index
    """
    return c * (r / r0)**gamma

# but to fit, may need to exclude bins where r < r0/2. or something like that to avoid non-linear scales. On large scale, it's set by primordial fluctuations that is ~ power-law.
maskii = np.argmin(np.abs(rbins - 1.0))

# Plot the Landy-Szalay TPCF and its fit
popt, pcov = curve_fit(TPCF_powerlaw, rbins[maskii:rsize-1],
                                      xi_all[maskii:rsize-1],
                        p0=[1., 2., 1.8])
plt.figure()
plt.plot(rbins[1:rsize-1], TPCF_powerlaw(rbin_centers[1:rsize-1], \
                                    popt[0], popt[1], popt[2]
                                    ), \
         ':k', linewidth=1,
         label=r'$\\xi(r)$=%.2f (r/%.2f)$^{%.2f}$'%(popt[0], popt[1], popt[2])
        )
plt.plot(rbin_centers[1:rsize-1], xi_all[1:rsize-1], label='All galaxies',
         color='k', marker='o')

plt.xlabel(r'$r $  $\rm{[Mpc h$^{-1}$]}$', fontsize=25)
plt.ylabel(r'$\\xi(r)$', fontsize=25)
plt.loglog()
plt.legend(loc='best', fontsize=20)
plt.savefig('realspace_xi_fitted_' + nn + '.pdf', bbox_inches='tight')


def construct_random_cat(Ndata, boxsize):

#    x = np.random.random(size=int(Ndata))
#    y = np.random.random(size=int(Ndata))
#    z = np.random.random(size=int(Ndata))

    x = np.random.uniform(0, boxsize, size=int(Ndata)*3)
    y = np.random.uniform(0, boxsize, size=int(Ndata)*3)
    z = np.random.uniform(0, boxsize, size=int(Ndata)*3)
    pos = np.vstack((x, y, z)).T
    return pos


r1 = construct_random_cat(len(galpos), boxsize)
xi, xi_cov = tpcf_jackknife(galpos, randoms=r1, rbins=rbins,
                            period=boxsize, Nsub=3)   # Nsub=5, Nsub=8

pfit_curvefit = xi
perr_curvefit = np.sqrt(np.absolute(np.diagonal(xi_cov)))     # np.array(error)


plt.figure()
plt.plot(rbin_centers, xi,
         label='All galaxies', color='k')
plt.errorbar(rbin_centers, pfit_curvefit, perr_curvefit, \
             fmt='.k', ecolor='gray', lw=1, label='')

plt.xlim(xmin=0.1, xmax=10)
plt.ylim(ymin=1, ymax=1e3)
plt.loglog()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r'$r $  $\rm{[Mpc]}$', fontsize=25)
plt.ylabel(r'$\\xi_{\rm gg}(r)$', fontsize=25)
plt.legend(loc='best', fontsize=20)
plt.savefig('halo_jackknife_' + nn + '.pdf', bbox_inches='tight')


# input array of host halo IDs that are equal for galaxies occupying the
# same halo
group_list = sim.galaxies
phid = []
for ii, o in enumerate(group_list):
    if o.halo is not None:
        phid.append(o.halo.GroupID)
    else:
        import pdb; pdb.set_trace()

print(len(galpos), len(phid))
xi_1h, xi_2h = tpcf_one_two_halo_decomp(galpos,
                                        phid, rbins,
                                        period=boxsize,
                                        num_threads='max')

plt.figure()
plt.plot(rbin_centers, xi_all,
         label='All galaxies', color='k')
plt.plot(rbin_centers, xi_1h,
         label='1-halo term')
plt.plot(rbin_centers, xi_2h,
         label='2-halo term')

plt.xlim(xmin=0.1, xmax=10)
plt.ylim(ymin=1, ymax=1e3)
plt.loglog()
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel(r'$r $  $\rm{[Mpc]}$', fontsize=25)
plt.ylabel(r'$\\xi_{\rm gg}(r)$', fontsize=25)
plt.legend(loc='best', fontsize=20)
plt.savefig('halo_decomp_' + nn + '.pdf', bbox_inches='tight')



# The integral of the TPCF as a function of cell volume --> cosmic variance

