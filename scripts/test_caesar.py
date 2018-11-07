'''

testing different funcs in caesar

'''


import caesar
import yt
import numpy as np
import pandas as pd

infile = '../simdata/s48_m25n256_z0.hdf5'
obj = caesar.load(infile)
redshiftFile = '../simdata/outputs_boxspace50.info'
_, zs_table, snaps_table = np.loadtxt(redshiftFile, unpack=True)

print obj.haloinfo(top=10)     # most massive halos
for i in range(100):
    print obj.halos[i].masses['total']
    # halos are sorted by total mass already

print len(obj.halos)
print len(obj.galaxies)        # total galaxies
print obj.ngalaxies

galaxy_masses = [i.masses['total'] for i in obj.galaxies]
print galaxy_masses
parent_halo_masses = [i.halo.masses['total'] for i in obj.galaxies]
big_galaxy_halo_masses = [i.halo.masses['total']
                          for i in obj.galaxies if i.masses['total'] > 1.0e12]
print big_galaxy_halo_masses
central_galaxy_halo_masses = [i.halo.masses['total']
                              for i in obj.galaxies if i.central]
Ngal = (np.array(central_galaxy_halo_masses) > 1.0e12).sum()
print("Number of central galaxies with masses >1.e12 Msun: {}").format(Ngal)
print("Corresponding to {:.2f}% of all central galaxies").format(
    Ngal * 100. / len(central_galaxy_halo_masses))

# co-moving, with unit cm attached
print obj.galaxies[0].radii['total']
print obj.galaxies[0].radii['total'].to('kpc')    # phyiscal kpc
print obj.galaxies[0].radii['total'].to('kpc/h')  # physical kpc/h
