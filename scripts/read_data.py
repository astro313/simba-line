'''
load in example sim data to see what the columns are

'''

import pandas as pd
# the z0 particle sim data only exist in release-branch
try:
    simdm = pd.read_pickle("temp/z0_data_files/particle_data/sim_data/z0.00_h401_s135_G2_sim.dm")
    simstar = pd.read_pickle("temp/z0_data_files/particle_data/sim_data/z0.00_h401_s135_G2_sim.star")
    simgas = pd.read_pickle("temp/z0_data_files/particle_data/sim_data/z0.00_h401_s135_G2_sim.gas")
    # The gas data must include: x, y, z, vx, vy, vx, mass, density, temperature, Z (at least an overall metallicity), molecular/dense gas mass fraction and instantaneous star formation rate.
except:
    print("These files only existing in release-branch of SIGAME_Dev")

print simdm.columns
# Index([u'm', u'vx', u'vy', u'vz', u'x', u'y', u'z'], dtype='object')

print simstar.columns
# Index([u'Z', u'age', u'm', u'vx', u'vy', u'vz', u'x', u'y', u'z', u'type',
      #  u'L_FUV'],
      # dtype='object')

print simgas.columns
# Index([u'SFR', u'Tk', u'Z', u'a_C', u'a_Ca', u'a_Fe', u'a_He', u'a_Mg', u'a_N',
      #  u'a_Ne', u'a_O', u'a_S', u'a_Si', u'f_H2', u'f_H21', u'f_HI', u'f_HI1',
      #  u'f_neu', u'h', u'm', u'nH', u'vx', u'vy', u'vz', u'x', u'y', u'z',
      #  u'type', u'FUV', u'CR', u'P_ext', u'surf_gas', u'surf_star',
      #  u'sigma_gas', u'sigma_star', u'vel_disp_gas'],
      # dtype='object')


# run extract_galaxy() on hdf5 files


# NOTE, in reality, sim data needs to reside in SIGAME_dev/sigame/temp/zx_data_files/particle_data/sim_data/ for SIGAME to run properly

# and a pickle file in temp/galaxies/ like z0_extracted_gals


# preparing for CLOUDY grid for redshift, in master branch:
# To calculate line emission from a raw galaxy dataset, SÍGAME needs to go through a set of tasks:
# 1. subgrid
#   - adding FUV field and external cloud pressure to the data
#   - creating GMCs and diffuse gas clouds
# --> files stored in particle_data/ISM_data/
# 2. interpolate
#   - interpolate cloud modesl to sum up line emission from GMC, DNG, and DIG ISM phases
# 3. datacubes
#   - drizzle the GMC, DNG, DIG results onto datacube
change parameters_zx.txt

import sigame as si
si.plot.grid_parameters_checkParam()

edit and run dif_cloud_grid.py
edit and run GMC_cloud_grid.py

si.run()


# make mom0 map from cube
si.plot.map_line()

