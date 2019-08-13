
"""

Number counts

"""


import caesar
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

# define input file
infile = '/mnt/ceph/users/daisyleung/simba/sim/m50n1024/s50/Groups/m50n1024_036.hdf5'

# load in input file
sim = caesar.load(infile,LoadHalo=False)
mlim = 64*sim.simulation.critical_density.value*sim.simulation.boxsize.value**3*sim.simulation.omega_matter/sim.simulation.effective_resolution**3 # galaxy mass resolution limit in DM particle mass
redshift = sim.simulation.redshift

ids = np.asarray([i.GroupID for i in sim.galaxies])
ms = np.asarray([i.masses['stellar'] for i in sim.galaxies])
mHI = np.asarray([i.masses['HI'] for i in sim.galaxies])
mdust = np.asarray([i.masses['dust'] for i in sim.galaxies])
mbh = np.asarray([i.masses['bh'] for i in sim.galaxies])
sfr = np.asarray([i.sfr for i in sim.galaxies])
met = np.asarray([i.metallicities['sfr_weighted'] for i in sim.galaxies])

S850 = 0.81 * ((sfr+1.e-6)/100)**0.43 * (mdust/1.e8)**0.54  # Hayward+11 fit to iso sims+Sunrise

fig,ax = plt.subplots()

n, bins, patches = ax.hist(S850, 20,
                            # normed=True,
                            histtype='step',
                            cumulative=-1,      # such that faintest bin would be sum of all points
                            label='Empirical')
ax.set_yscale('log')

plt.minorticks_on()
plt.ylabel(r'$\log\ N(>S)$',fontsize=16)
plt.xlabel(r'$\log\ S_{850}$' ,fontsize=16)
plt.annotate('z=%g'%(np.round(redshift,1)), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))

plt.show()

