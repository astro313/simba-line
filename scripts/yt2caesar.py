"""

Regardless of RD's comment on caesar gas field, we are sticking with
caesar gas mass so that gas_f_neu, f_h2 is consistent.
Otherwise the number don't add up. Meaning that

gal.masses['HI'] >> np.sum(gas_f_neu * gas_m)
so then gas_f_neu wouldn't be useful anymore.

"""


from __future__ import print_function, division
from astropy import constants as constants
from readgadget import *

# for consistency, will use python 3 for all scripts of this project.
import sys
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

import _pickle as pickle
from parse_simba import *

import os
import sys
import numpy as np
import caesar
import matplotlib.pyplot as plt


def get_partmasses_from_snapshot(snapFile, obj, ptype, physicalUnit=True, verbose=False):
    """

    Use readgadget to read in the info of each particles.

    Parameters
    ----------
    snapFile: str

    obj: caesar obj

    ptype: str
        'dm', 'star', 'gas', 'dust'

    physicalUnits: bool
        should be True in most cases

    Returns
    -------
    m: array
        array of mass for each particle

    NOTE
    ----
    Note "gas" mass is just the "cold gas" identified in the galaxy using the 6dFOF group finding which is restricted to cold gas > 0.13 atom/cc, does NOT include all the HI since self-shielding kicks in at ~1E-3 cc.

    RD: usually ignore 'gas' and just directly use 'HI' and 'H2' as total gas content.


    """

    if physicalUnit:
        units = 1
    else:
        units = 0

    h = obj.simulation.hubble_constant

    if ptype == 'dm':
        try:
            m = readsnap(snapFile,'mass','dm',units=units)/h
        except:
            return 0.0

    elif ptype == 'gas':
        m = readsnap(snapFile,'mass','gas',units=units)/h    # solar mass

    elif ptype == 'dust':
        m = readsnap(snapFile,'Dust_Masses','gas')*1.e10/h

    elif ptype == 'star':
        m = readsnap(snapFile,'mass','star',units=units,
                          suppress=1)/h  # suppress means no output print to command line
    if verbose:
        print("Number of {:} particles in snapshot file {:}: {:}").format(ptype, snapFile, len(m))

    return m


def group_part_by_galaxy(snapPart, galaxy, ptype):
    """
        take an array of particles and then return an array containing only the particles that belong to the given galaxy.
    """

    if ptype == 'gas':
        snapPart = np.array([snapPart[k] for k in galaxy.glist])
    elif ptype == 'star':
        snapPart = np.array([snapPart[k] for k in galaxy.slist])
    else:
        raise NotImplemented("Unclear ptype")
    return snapPart



def check_dense_gas(dir='./'):
    """

    Check the number of gas particles in each galaxy that has non-zero H2 fraction. This will determine how many GMCs we will get in the subgrid step in galaxy.add_GMCs(). If this number is low, then we might as well skip extracting that galaxy.

    This should no longer be needed after we've added the denseGasThres in the particles2pd class.

    """
    import glob
    import pandas as pd
    ff = glob.glob('*.gas')

    for i in ff:
        f = pd.read_pickle(i)
        print(i)
        print (f['f_H21'] > 0.0).sum()

        print("Total dense gas mass: ")
        print(f['m'] * f['f_H21']).sum()
    return None


class particles2pd(object):


    def __init__(self, snapRange=[36], name_prefix='m25n1024_', feedback='s50/', zCloudy=6, part_threshold=64, sfr_threshold=0.1, denseGasThres=1.e5, user='Daisy', debug=False, verbose=True):
        """

        Parameters
        ----------
        snapRange: list of int
            list of snapshots to use

        name_prefix: str
            determines which simulation box to read in

        feedback: str
            determines which simulation box to read in

        zCloudy: float
            redshift to use for cloudy table

        part_threshold: int
            both stars and gas particles must be > part_threshold particles each.

        sfr_threshold: float
            lower limit of SFR to include in sample

        denseGasThres: float
            min. gas in dense phase in Msun in order for a galaxy to be included. If too small, we won't be able to do GMC subgridding.
        """

        self.Mp = 1.67262189821e-24
        self.kpc2m = 3.085677580666e19

        self.snapRange = snapRange    # don't put 036
        self.zCloudy = zCloudy
        self.name_prefix = name_prefix

        if feedback[-1] != '/': feedback += '/'
        self.feedback = feedback

        self.raw_sim_name_prefix = 'snap_' + name_prefix
        self.part_threshold = part_threshold
        self.sfr_threshold = sfr_threshold
        self.denseGasThres = denseGasThres       # Msun

        self.user = user

        self.debug = debug
        self.verbose = verbose
        self.setup()

    def setup(self):
        """

        setup path to simulation file and path to save pandas dataframe.
        Only implemented for Daisy, works on ursa or the Flatiron Server.

        """

        if self.user is 'Daisy':
            import socket
            host = socket.gethostname()

            simName = self.name_prefix[:self.name_prefix.find('_')]

            if 'ursa' in host:
                self.raw_sim_dir = '/disk01/rad/sim/' + simName + '/' + self.feedback
                self.caesar_dir = '/disk01/rad/sim/' + simName + '/' + self.feedback + 'Groups/'
                self.redshiftFile = '/home/rad/gizmo-extra/outputs_boxspace50.info'
                self.d_data = '/home/dleung/Downloads/SIGAME_dev/sigame/temp/z' + str(int(self.zCloudy)) + '_data_files/'
            elif 'flatironinstitute.org' or 'worker' in host:
                self.raw_sim_dir = '/mnt/ceph/users/daisyleung/simba/sim/' + simName + '/' + self.feedback  # dummy
                self.caesar_dir = '/mnt/ceph/users/daisyleung/simba/sim/' + simName + '/' + self.feedback + 'Groups/'
                self.redshiftFile = '/mnt/ceph/users/daisyleung/simba/gizmo-extra/outputs_boxspace50.info'
                self.d_data = '/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/z' + str(int(self.zCloudy)) + '_data_files/'
        else:
            raise NotImplementedError


    def def_caesarFileName(self):
        infile = self.caesar_dir + self.name_prefix + '{:0>3}'.format(int(self.snap)) + \
            '.hdf5'
        return infile

    def def_snapFileName(self):

        snapfile = self.raw_sim_dir + self.raw_sim_name_prefix + '{:0>3}'.format(int(self.snap)) + '.hdf5'
        return snapfile


    def load_obj_snap(self, idx, redshiftDecimal=2, LoadHalo=False):
        self.snap = self.snapRange[idx]

        infile = self.def_caesarFileName()
        print("Loading Ceasar file: {:}".format(infile))
        self.obj = caesar.load(infile, LoadHalo=LoadHalo)

        self.h = self.obj.simulation.hubble_constant
        # self.redshift = np.round(self.obj.simulation.redshift, redshiftDecimal)
        self.redshift = self.obj.simulation.redshift

        # snap
        self.snapFile = self.def_snapFileName()


    def run(self, savepath=None, outname=None, emptyDM=True, caesarRotate=False):
        """
        Loop through snapRange and run main_proc()

        Parameters
        ----------
        savepath: str
            where the pandas .gas, .star, .dm should be saved to
            default is None, which means subdir of sigame

        emptyDM: bool
            if True, will create empty holders as DM position, mass, velocity to save as pandas dataframe. Do so because we don't have a dmlist from caesar catalog and we don't need it really.

        outname: str
            path to where zx_extracted_gals is stored for global_results.py to read

        caesarRotate: bool
            whether to project gal to xy-plane.

        """

        for idx in range(len(self.snapRange)):
            self.load_obj_snap(idx)
            if idx == 0:
                gnames, zzz = self.main_proc(savepath=savepath, emptyDM=emptyDM, caesarRotate=caesarRotate)
            else:
                gnamesOut, zzzRed = self.main_proc(savepath=savepath, emptyDM=emptyDM, caesarRotate=caesarRotate)
                gnames.extend(gnamesOut)
                zzz.extend(zzzRed)

        from parse_simba import pd_bookkeeping
        _, _ = pd_bookkeeping(ganmes, zzz, self.zCloudy, outname=outname)

        return gnames, zzz


    def main_proc(self, savepath, emptyDM, caesarRotate):

        """
        Parameters
        ----------

        savepath: str
            where the pandas .gas, .star, .dm should be saved to
            default is None, which means subdir of sigame

        emptyDM: bool
            if True, will create empty holders as DM position, mass, velocity to save as pandas dataframe. Do so because we don't have a dmlist from caesar catalog and we don't need it really.

        caesarRotate: bool
            whether to project gal to xy-plane.

        Returns
        -------
        galName: list of str
            galnames extracted, after applying the selection criteria

        Note
        ----
        We will apply some selection criteria:
        - SFR
        - dense gas mass, otherwise we will get an error in subgrid add_GMCs()
        - stars and gas particles number


        """

        from yt2caesar import get_partmasses_from_snapshot
        from yt2caesar import group_part_by_galaxy
        from astropy import units as u
        import pandas as pd

        if savepath is None:
            savepath = self.d_data + 'particle_data/sim_data/'

        # sort by SFR
        self.obj.galaxies.sort(key=lambda x: x.sfr, reverse=True)

        # load in the fields from snapshot
        print("Read in gas fields")

        if self.debug:
            # compare against output saved from YT (parse_simba.py)
            sim_gas = pd.read_pickle('/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/z6_data_files/particle_data/sim_data/z5.93_h0_s36_G0_sim.gas')
            sim_star = pd.read_pickle('/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/z6_data_files/particle_data/sim_data/z5.93_h0_s36_G0_sim.star')

        rho_crit_cgs = 1.8791e-29      # /h^2
        unit_Density = rho_crit_cgs *self.h*self.h * u.g/(u.cm**3)
        gas_densities_p = readsnap(self.snapFile,'rho','gas',units=1)
                          # gas density in comoving g/cm^3
        # print("density g/cc")
        # print(gas_densities_p.min(), gas_densities_p.max())
        #

        if self.verbose or self.debug:
            # not acutally used in the dataframe
            gas_nh_p = gas_densities_p*self.h*self.h*0.76/self.Mp     # number density of H in 1/cc

        gas_p_m = get_partmasses_from_snapshot(self.snapFile, self.obj,
                                               ptype='gas')
        gas_Tk_p = readsnap(self.snapFile,'u','gas',units=1)
        gas_SFR_p = readsnap(self.snapFile,'sfr','gas',units=1, \
                             suppress=1)/self.h
        gmet_p = readsnap(self.snapFile,'metals','gas',units=1)
        # Smoothing length
        gas_h_p = readsnap(self.snapFile,'hsml','gas',units=1) \
                           /self.h/(1+self.redshift)  # smoothing length in ckpc --> proper kpc
        #
        gas_pos_p = readsnap(self.snapFile,'pos','gas',units=1,suppress=1)\
                             /self.h # ckpc
        gas_vel_p = readsnap(self.snapFile,'vel','gas',units=1,suppress=1)   # km/s

        # molecular gas fraction
        gfH2_p = readsnap(self.snapFile,'fH2','gas', units=1)
        assert abs(gfH2_p.all()) <= 1.0      # each particles
        # neutral hydrogen fraction (between 0-1)
        gfHI_p = readsnap(self.snapFile,'nh','gas',units=1)
        assert abs(gfHI_p.all()) <= 1.0      # each particles

        gas_x_e_p = readsnap(self.snapFile,'ne','gas',units=1)

        if self.debug:
            nH = gfHI_p*gas_densities_p/self.Mp              # number density
            print('nH: ')
            print('from gfHI and gas density of readsnap')
            print(nH.max(), nH.min())
            print("from gas density of readsnap and 0.76: ")
            print(gas_nh_p.max(), gas_nh_p.min())
            plt.figure()
            plt.hist(nH, bins=100)
            plt.hist(gas_nh_p, bins=100)
            plt.show(block=False)
            plt.savefig('123.pdf')
            nE = gas_x_e_p * nH
           # nE = [a*b for a,b in zip(gas_x_e_p, nH)]
            print("Electron Number Density: ")
            print(nE)
            import pdb; pdb.set_trace()

        print("Read in stellar fields")
        star_p_m = get_partmasses_from_snapshot(self.snapFile, self.obj, ptype='star')
        star_pos_p =readsnap(self.snapFile,'pos','star',units=1,suppress=1)/self.h   # ckpc
        star_vel_p =readsnap(self.snapFile,'vel','star',units=1,suppress=1)
        pmetarray = readsnap(self.snapFile,'Metallicity','star',units=1,suppress=1)[:, 0]
        sage = readsnap(self.snapFile,'age','star',units=1,suppress=1)    #  expansion factor of formation

        print("Read in DM mass for each particle: " )
        dm_p_m = get_partmasses_from_snapshot(self.snapFile, self.obj, ptype='dm')

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        SolarAbundances=[0.0134, 0.2485, 2.38e-3, 0.70e-3, 5.79e-3, 1.26e-3,
                         7.14e-4, 6.17e-4, 3.12e-4, 0.65e-4, 1.31e-3]

        galName = []
        zred = []
        for gg, gal in enumerate(self.obj.galaxies):

            loc = gal.pos    # ckpc

            galname = 'h' + str(int(gal.parent_halo_index)) + '_s' + \
                str(int(self.snap)) + '_G' + str(int(gg))

            if self.verbose:
                print(galname)
                print("SFR: {:.2f}".format(gal.sfr))

            gas_m = group_part_by_galaxy(gas_p_m, gal, ptype='gas')
            gas_densities = group_part_by_galaxy(gas_densities_p, gal, ptype='gas')

            if self.debug:
                print("from readsnap: ")
                print(gas_densities.max(), gas_densities.min())    # g/cc
                print("from YT sphere: ")
                print(sim_gas['nH'].max(), sim_gas['nH'].min())
                gas_nh = group_part_by_galaxy(gas_nh_p, gal, ptype='gas')
                print(gas_nh.max(), gas_nh.min())      # 1/cc
                import pdb; pdb.set_trace()

            gas_Tk = group_part_by_galaxy(gas_Tk_p, gal, ptype='gas')
            gas_SFR = group_part_by_galaxy(gas_SFR_p, gal, ptype='gas')
            gas_Z = group_part_by_galaxy(gmet_p[:, 0], gal, ptype='gas')/SolarAbundances[0]
            gas_Z_1 = group_part_by_galaxy(gmet_p[:, 1], gal, ptype='gas')/SolarAbundances[1]
            gas_Z_2 = group_part_by_galaxy(gmet_p[:, 2], gal, ptype='gas')/SolarAbundances[2]
            gas_Z_3 = group_part_by_galaxy(gmet_p[:, 3], gal, ptype='gas')/SolarAbundances[3]
            gas_Z_4 = group_part_by_galaxy(gmet_p[:, 4], gal, ptype='gas')/SolarAbundances[4]
            gas_Z_5 = group_part_by_galaxy(gmet_p[:, 5], gal, ptype='gas')/SolarAbundances[5]
            gas_Z_6 = group_part_by_galaxy(gmet_p[:, 6], gal, ptype='gas')/SolarAbundances[6]
            gas_Z_7 = group_part_by_galaxy(gmet_p[:, 7], gal, ptype='gas')/SolarAbundances[7]
            gas_Z_8 = group_part_by_galaxy(gmet_p[:, 8], gal, ptype='gas')/SolarAbundances[8]
            gas_Z_9 = group_part_by_galaxy(gmet_p[:, 9], gal, ptype='gas')/SolarAbundances[9]
            gas_Z_10 = group_part_by_galaxy(gmet_p[:, 10], gal, ptype='gas')/SolarAbundances[10]

            # smoothing length
            gas_h = group_part_by_galaxy(gas_h_p, gal, ptype='gas')

            gas_pos = group_part_by_galaxy(gas_pos_p, gal, ptype='gas')

            if len(gal.slist) < self.part_threshold or len(gal.glist) < self.part_threshold:
                print("Too few star particles or gas particles, unlikely to be real galaxy or useful for our purpose. Skipping ", galname)
                continue

            gas_pos -= loc.d          # both are in comoving
            gas_pos /= (1+self.redshift)   # physical kpc
            gas_vel = group_part_by_galaxy(gas_vel_p, gal, ptype='gas')

            if caesarRotate:
                gas_pos = caesar.utils.rotator(gas_pos, gal.rotation_angles['ALPHA'].astype('float64'), gal.rotation_angles['BETA'].astype('float64'))
                gas_vel = caesar.utils.rotator(gas_vel,
                                              np.float64(gal.rotation_angles['ALPHA']),
                                              np.float64(gal.rotation_angles['BETA']))
                ff = lambda x: x.d
                gas_x = ff(gas_pos[:, 0])
                gas_y = ff(gas_pos[:, 1])
                gas_z = ff(gas_pos[:, 2])

                gas_vx = ff(gas_vel[:,0])
                gas_vy = ff(gas_vel[:,1])
                gas_vz = ff(gas_vel[:,2])
            else:
                gas_x = gas_pos[:, 0]
                gas_y = gas_pos[:, 1]
                gas_z = gas_pos[:, 2]

                gas_vx = gas_vel[:,0]
                gas_vy = gas_vel[:,1]
                gas_vz = gas_vel[:,2]

            gas_f_H2 = group_part_by_galaxy(gfH2_p, gal, ptype='gas')
            gas_f_neu = group_part_by_galaxy(gfHI_p, gal, ptype='gas')   # f_neu

            if self.debug:
                print("HI fraction")
                print(gas_f_neu.min(), gas_f_neu.max())    # following RD's def.
                print((1-gas_f_H2).min(), (1-gas_f_H2).max())
                print((gas_f_neu/(1-gas_f_H2)).min())
                print((gas_f_neu/(1-gas_f_H2)).max())

            # neutral gas from 1- ionized gas
            gas_x_e = group_part_by_galaxy(gas_x_e_p, gal, ptype='gas')   # relative to nH
            gas_f_ion = gas_x_e / max(gas_x_e)
            gas_f_HI = 1 - gas_f_ion
            assert abs(gas_f_HI.all()) <= 1.0

            if self.debug:
                print(gas_f_HI >= gas_f_neu)     # former incl. also molecular
                print((gas_f_HI >= gas_f_neu).all())   # expecting True
                import pdb; pdb.set_trace()

            if self.debug:
                print('\nChecking molecular gas mass fraction from simulation:')
                print('%.3s %% \n' % (np.sum(gas_m * gas_f_H2) / np.sum(gas_m) * 100.))
                #
                print("gas mass from snapshot: {:.2f} [x1e8 Msun]".format(gas_m.sum()/1.e8))
                print("gas mass from 'gas' from caesar {:.2f} [x1e8 Msun]".format(gal.masses['gas']/1.e8))
                #
                print('gas mass from (HI + H2) from caesar {:.2f} [x1e8 Msun]'.format((gal.masses['HI'] + gal.masses['H2'])/1.e8))
                print('')
                print("gas fraction from caesar: {:.2f}".format(gal.gas_fraction))
                print('gas fraction from Mgas/(Mgas+Mstar): {:.2f} '.format(gal.masses['gas']/(gal.masses['gas'] + gal.masses['stellar'])))
                print('gas fraction from MHI + MH2 /(MHI + MH2 + Mstar): {:.2f}'.format((gal.masses['HI'] + gal.masses['H2']) / (gal.masses['HI'] + gal.masses['H2'] + gal.masses['stellar'])))
                #
                print(gal.masses['HI'], np.sum(gas_f_neu * gas_m))
                import pdb; pdb.set_trace()

            # selection crit.
            if gas_SFR.sum() <= self.sfr_threshold:
                print("SFR too low.. Skipping ", galname)
                continue

            if (gas_m * gas_f_H2).sum() <= self.denseGasThres:
                print ("Dense gas mass less than %.2f Msun.. Skipping %s" % (denseGasThres, galname ))
                continue

            if len(gal.slist) < self.part_threshold or len(gal.glist) < self.part_threshold:
                print("Too few star particles or gas particles, unlikely to be real galaxy or useful for our purpose. Skipping ", galname)
                continue

            star_m = group_part_by_galaxy(star_p_m, gal, ptype='star')
            star_pos = group_part_by_galaxy(star_pos_p, gal, ptype='star')
            star_pos -= loc.d
            star_vel = group_part_by_galaxy(star_vel_p, gal, ptype='star')

            if caesarRotate:
                star_pos = caesar.utils.rotator(star_pos,
                                   np.float64(gal.rotation_angles['ALPHA']),
                                   np.float64(gal.rotation_angles['BETA']))

                star_vel = caesar.utils.rotator(star_vel,
                                    np.float64(gal.rotation_angles['ALPHA']),
                                    np.float64(gal.rotation_angles['BETA']))

                star_x = ff(star_pos[:, 0])
                star_y = ff(star_pos[:, 1])
                star_z = ff(star_pos[:, 2])

                star_vx = ff(star_vel[:,0].d)
                star_vy = ff(star_vel[:,1].d)
                star_vz = ff(star_vel[:,2].d)
            else:
                star_x = star_pos[:, 0]
                star_y = star_pos[:, 1]
                star_z = star_pos[:, 2]

                star_vx = star_vel[:,0]
                star_vy = star_vel[:,1]
                star_vz = star_vel[:,2]

            star_Z =  group_part_by_galaxy(pmetarray, gal, ptype='star')

            # derive stellar age
            star_a = group_part_by_galaxy(sage, gal, ptype='star')
            current_time = self.obj.simulation.time.in_units("Myr")
            # in scale factors, do as with Illustris
            star_formation_z = 1. / star_a - 1
            # Code from yt project (yt.utilities.cosmology)
            star_formation_t = 2.0 / 3.0 / np.sqrt(1 - self.obj.simulation.omega_matter) * np.arcsinh(np.sqrt(
                (1 - self.obj.simulation.omega_matter) / self.obj.simulation.omega_matter) / np.power(1 + star_formation_z, 1.5)) / (self.h)  # Mpc*s/(100*km)
            star_formation_t *=  self.kpc2m / 100. / (1e6 * 365.25 * 86400)  # Myr
            star_age = current_time.d - star_formation_t

            # create empty DM data
            if emptyDM:
                dm_m = 0.0

                dm_posx = np.array([0.0])
                dm_posy = np.array([0.0])
                dm_posz = np.array([0.0])

                dm_velx = np.array([0.0])
                dm_vely = np.array([0.0])
                dm_velz = np.array([0.0])
            else:
                # because caesar output dones't have a dmlist and since we don't really need it
                raise NotImplementedError

            # create pandas DF
            simgas_path = (savepath + 'z{:.2f}').format(float(self.redshift)) + '_' + \
                           galname + '_sim.gas'
            simstar_path = (savepath + 'z{:.2f}').format(float(self.redshift)) + \
                           '_' + galname + '_sim.star'
            simdm_path = (savepath + 'z{:.2f}').format(float(self.redshift)) + \
                          '_' + galname + '_sim.dm'

            simgas = pd.DataFrame({'x': gas_x, 'y': gas_y, 'z': gas_z,
                                   'vx': gas_vx, 'vy': gas_vy, 'vz': gas_vz,
                                   'SFR': gas_SFR, 'Z': gas_Z,
                                   'nH': gas_densities,
                                   'Tk': gas_Tk, 'h': gas_h,
                                   'f_HI1': gas_f_HI,    # atomic and molecular H
                                   'f_neu': gas_f_neu,   # atomic H
                                   'f_H21': gas_f_H2, 'm': gas_m,
                                   'a_He': gas_Z_1, 'a_C': gas_Z_2,
                                   'a_N': gas_Z_3, 'a_O': gas_Z_4,
                                   'a_Ne': gas_Z_5, 'a_Mg': gas_Z_6,
                                    'a_Si': gas_Z_7, 'a_S': gas_Z_8,
                                    'a_Ca': gas_Z_9, 'a_Fe': gas_Z_10})

            simstar = pd.DataFrame({'x': star_x, 'y': star_y, 'z': star_z,
                                    'vx': star_vx, 'vy': star_vy, 'vz': star_vz,
                                    'Z': star_Z, 'm': star_m, 'age': star_age})

            # create fake DM dataframe to trick Sigame
            simdm = pd.DataFrame({'x': dm_posx, 'y': dm_posy, 'z': dm_posz,
                                  'vx': dm_velx, 'vy': dm_vely, 'vz': dm_velz,
                                  'm': dm_m})


            # from parse_simba import center_cut_galaxy
            # simgas, simstar, simdm = center_cut_galaxy(simgas, simstar, simdm, plot=False)
            # import pdb; pdb.set_trace()


            simgas.to_pickle(simgas_path)
            simstar.to_pickle(simstar_path)
            simdm.to_pickle(simdm_path)

            galName.append(galname)
            zred.append(self.redshift)
        return galName, zred



if __name__ == '__main__':

    pp = particles2pd(snapRange=[36],name_prefix='m25n1024_', feedback='s50/', zCloudy=6, user='Daisy', part_threshold=64, sfr_threshold=0.1, denseGasThres=1.e4)

    ggg, zred = pp.run(savepath='xxx/', outname=None, emptyDM=True, caesarRotate=False)
    print(ggg)



''' Ways to select physically meaningful galaxies ....

ghaloall = np.array(readsnap(snap,'HaloID','gas',suppress=1,nth=nskip),dtype=int)  # FOF halo ID from snapshot
gas_select = ((gnh>0.1)&(ghaloall>0))  # selects gas to use

if gal.masses['bh'] <= 1.e6: continue    # a way to select gal w/ BH


'''
