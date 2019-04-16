'''

parse halo-galaxy catalog from caesar output file and get what i need.

'''

def get_basic_info_from_caesarCat(snapRange, Nhalo, Ngalaxies, caesar_dir, name_prefix):

    """

    Read the caesar file and get some basic info on num of galaxy in those snapshots

    """

    import caesar
    for ii, sss in enumerate(snapRange):
        infile = caesar_dir + name_prefix + '{:0>3}'.format(int(sss)) + \
            '.hdf5'
        print("Loading Ceasar file: {}").format(infile)
        obj = caesar.load(infile)

    print('Total number of galaxies found: ' + str(obj.ngalaxies))
    Ngal = obj.ngalaxies

    print('Info on the massive {} halos: '.format(Nhalo))
    print(obj.haloinfo(top=Nhalo))

    print('Info on the massive {} galaxies across all haloes: '.format(Ngalaxies))
    print(obj.galinfo(top=Ngalaxies))

    central_galaxy_halo_masses = [i.halo.masses['total']
                                  for i in obj.galaxies[:Ngalaxies] if i.central]
    print('Where {}% are centrals'.format(len(central_galaxy_halo_masses) * 100. / Ngalaxies))

    print("The top Ngalaxies most massive galaxies reside in halos ID: ")
    for ggg in range(Ngalaxies):
        # all fields:
        # print(obj.galaxies[ggg].info())
        print("Halo ID: {}".format(obj.galaxies[ggg].parent_halo_index))

    print("\nCount the number of galaxies in each massive halo among top Nhalo: ")
    for hhh in range(Nhalo):
        numGal = len(obj.halos[hhh].galaxies)
        print("Found total of {} galaxies in halo {}".format(numGal, hhh))

    return Ngal


def select_SFgal_from_simba(raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, snapRange, Ngalaxies, saveggg=None, verbose=False, debug=False):

    '''
        pick out the 'Ngalaxies' most star-forming galaxies across snapshots 'snapRange'

    Parameters
    ----------
    raw_sim_dir: string
        where to look for raw snapshot files
    raw_sim_name_prefix: string
        prefix to the speciifc volume and resolution of snapshots we are looking
    caesar_dir: str
        path to the caesar output files
    name_prefix: str
        predix to file of caesar outputs
    snapRange: list of int
        snapshots to look for galaxies
    Ngalaxies: int
        how many galaxies from each snapshot across all halos do we want as output
    saveggg: str
        filename to save output "galnames_selected"

    Returns
    -------
    galnames_selected: list of string
        galaxies names indicating the halo ID, snapshot number, and galaxy ID (numbered based on some predefined criterion)

    '''

    import caesar
    import numpy as np
    import pandas as pd

    # store gal info for large number of galaxies, which we then select from
    galnames = pd.DataFrame({'halo': [], 'snap': [], 'GAL': [], 'SFR': []})

    # loop thru snapshots
    for ii, sss in enumerate(snapRange):
        infile = caesar_dir + name_prefix + '{:0>3}'.format(int(sss)) + \
            '.hdf5'
        print("Loading Ceasar file: {}").format(infile)
        obj = caesar.load(infile)    # # LoadHalo=False

        print("\nSelecing {} Galaxies with the highest SFRs across all halos in this snapshot, but you may want to galaxies based on different criteria.").format(Ngalaxies)
        obj.galaxies.sort(key=lambda x: x.sfr, reverse=True)
        if verbose:
            print(obj.galinfo(top=Ngalaxies))

        for GAL in range(Ngalaxies):
            try:
                SFR = float(obj.galaxies[GAL].sfr.d)
                add_this = pd.DataFrame({'halo': [int(obj.galaxies[GAL].parent_halo_index)], 'snap': [
                                        sss], 'GAL': [GAL], 'SFR': [SFR]})
                galnames = galnames.append(add_this, ignore_index=True)
                galnames[['halo', 'snap', 'GAL']] = galnames[
                    ['halo', 'snap', 'GAL']].astype(int)
            except:
                break

    print("\nSelecing {} Galaxies with the highest SFRs across snapshot 'snapRange', but you may want to galaxies based on different criteria.\nIf so, edit define_ds()").format(Ngalaxies)
    galnames = galnames.sort_values(['SFR'], ascending=False).reset_index(drop=True)
    if debug:
        print(galnames)
    print("Note to self: Remember to change variable global_save_files in param.py so that naming is consistent with Ngalaxies we are picking here")
    galnames = galnames[:Ngalaxies]
    print(galnames)

    if saveggg is not None:
        with open(saveggg, 'w') as f:
            pickle.dump([ggg], f)
    return galnames


def print_ds_info(h, obj, ds):
    boxsize = (obj.simulation.boxsize*h).d       # boxsize=readheader(snap,'boxsize')
    # galaxy mass resolution limit: 32 gas particle masses
    mlim = 32*obj.simulation.critical_density.value*obj.simulation.boxsize.value**3*obj.simulation.omega_baryon/obj.simulation.effective_resolution**3
    print("Box size: {:} Mpc/h").format(boxsize/1.e3)
    print("Galaxy mass lower bound/mass reolution limit (32 particles masses): {:} ").format(mlim)
    print('Info for this snapshot:')
    print('Simulation type: ' + ds.dataset_type)
    for key in ds.parameters.keys():
        print('%s = %s' % (key, ds.parameters[key]))


def define_ds(snap, raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile, verbose=False):

    """

    """

    import caesar
    import yt
    import numpy as np
    import os
    import pygad as pg

    _, zs_table, snaps_table = np.loadtxt(redshiftFile, unpack=True)

    infile = caesar_dir + name_prefix + '{:0>3}'.format(int(snap)) + \
        '.hdf5'
    print("Loading Ceasar file: {}").format(infile)
    obj = caesar.load(infile)
    obj.galaxies.sort(key=lambda x: x.sfr, reverse=True)

    rawSim = raw_sim_dir + raw_sim_name_prefix + \
        '{:>03}'.format(int(snap)) + '.hdf5'

    ds = yt.load(rawSim, over_refine_factor=1, index_ptype="all")

    s = pg.Snap(raw_sim_dir + raw_sim_name_prefix + '{:0>3}'.format(int(snap)) + '.hdf5')
    h = s.cosmology.h()    # obj.simulation.hubble_constant

    if verbose:
        print_ds_info(h, obj, ds)

    zred = '{:.3f}'.format(zs_table[snaps_table == snap][0])
    return s, h, obj, ds, float(zred)


def center_cut_galaxy(simgas, simstar, simdm, plot=False):
    print(' ** Center galaxy in spatial and velocity coordinates, gas+stars, and cut out! **')
    plt.close('all')

    print('Number of gas particles: ' + str(len(simgas)))
    print('Number of star particles: ' + str(len(simstar)))
    print('Number of DM particles: ' + str(len(simdm)))

    # Define R_gal from simulation cutout
    R_gal = max(np.sqrt(simgas['x']**2. + simgas['y']**2. + simgas['z']**2.))
    print('R_gal = ' + str.format("{0:.2f}", R_gal) + ' kpc (from simulation)')

    print('Center all in x,y,z spatial coordinates (according to stellar distribution)')
    dat = simstar
    bins = 200.
    r = np.sqrt(dat['x'].values**2. + dat['y'].values**2. + dat['z'].values**2.)
    # r_bin           =   np.arange(-max(r),max(r),2.*max(r)/bins)
    r_bin = np.arange(-2, 2, 4. / bins)
    # Average mass surface density in radial bins:
    m_binned_x = np.array([sum(dat['m'][(dat['x'] >= r_bin[i]) & (
        dat['x'] < r_bin[i + 1])]) for i in range(0, len(r_bin) - 1)])
    m_binned_y = np.array([sum(dat['m'][(dat['y'] >= r_bin[i]) & (
        dat['y'] < r_bin[i + 1])]) for i in range(0, len(r_bin) - 1)])
    m_binned_z = np.array([sum(dat['m'][(dat['z'] >= r_bin[i]) & (
        dat['z'] < r_bin[i + 1])]) for i in range(0, len(r_bin) - 1)])
    # Smooth out profiles a bit:
    m_binned_x1 = lowess(m_binned_x, r_bin[0:len(
        r_bin) - 1], frac=0.1, is_sorted=True, it=0)
    m_binned_y1 = lowess(m_binned_y, r_bin[0:len(
        r_bin) - 1], frac=0.1, is_sorted=True, it=0)
    m_binned_z1 = lowess(m_binned_z, r_bin[0:len(
        r_bin) - 1], frac=0.1, is_sorted=True, it=0)
    # find max of distribution:
    xpos = r_bin[np.argmax(m_binned_x1[:, 1])]
    ypos = r_bin[np.argmax(m_binned_y1[:, 1])]
    zpos = r_bin[np.argmax(m_binned_z1[:, 1])]
    print('corrections: ', xpos, ypos, zpos)
    # move original coordinates
    simgas[pos] = simgas[pos] - [xpos, ypos, zpos]
    simstar[pos] = simstar[pos] - [xpos, ypos, zpos]
    simdm[pos] = simdm[pos] - [xpos, ypos, zpos]
    r_bin = r_bin[0:len(r_bin) - 1]
    if plot:
        plt.figure(0)
        plt.plot(r_bin, m_binned_x, 'r')
        plt.plot(r_bin, m_binned_y, 'g')
        plt.plot(r_bin, m_binned_z, 'b')
        plt.plot(r_bin, m_binned_x1[:, 1], '--r')
        plt.plot(r_bin, m_binned_y1[:, 1], '--g')
        plt.plot(r_bin, m_binned_z1[:, 1], '--b')
        # plt.plot([xpos,xpos],[0,1e10],'--r')
        # plt.plot([ypos,ypos],[0,1e10],'--g')
        # plt.plot([zpos,zpos],[0,1e10],'--b')
        plt.xlabel('x (r) y (g) z (b) [kpc]')
        plt.ylabel('accumulated stellar mass [M$_{\odot}$]')
        plt.title('Centering in [x,y,z]')
        plt.show(block=False)

    print('Center all in velocity space (vx,vy,vz) (according to gas distribution)')
    dat = simgas
    # use gas to center galaxy in velocity space (as if observing)
    ngrid = 1000
    grid, vxd, vyd, vzd = ([0] * ngrid for i in range(4))
    grid[0] = -600.
    for i in range(1, len(grid)):
        grid[i] = grid[i - 1] + 2 * (-grid[0]) / ngrid
        vxd[i] = sum(dat.loc[dat.loc[:, 'vx'] < grid[i], 'm'])
        vyd[i] = sum(dat.loc[dat.loc[:, 'vy'] < grid[i], 'm'])
        vzd[i] = sum(dat.loc[dat.loc[:, 'vz'] < grid[i], 'm'])
    # find the position where half of the mass has accumulated
    vxpos = max(np.array(grid)[np.array(vxd) < max(vxd) / 2])
    vypos = max(np.array(grid)[np.array(vyd) < max(vyd) / 2])
    vzpos = max(np.array(grid)[np.array(vzd) < max(vzd) / 2])
    # correct velocities
    simgas[vpos] = simgas[vpos] - [vxpos, vypos, vzpos]
    simstar[vpos] = simstar[vpos] - [vxpos, vypos, vzpos]
    # simdm[vpos]      =   simdm[vpos]-[vxpos,vypos,vzpos]
    print('corrections: ', vxpos, vypos, vzpos)


    r = rad(simgas[pos], pos)
    simgas = simgas[r < R_gal]
    simgas = simgas.reset_index(drop=True)
    r = rad(simstar[pos], pos)
    simstar = simstar[r < R_gal]
    simstar = simstar.reset_index(drop=True)
    r = rad(simdm[pos], pos)
    simdm = simdm[r < R_gal]
    simdm = simdm.reset_index(drop=True)

    # Make dataframe with global properties
    M_star = sum(simstar['m'])
    M_gas = sum(simgas['m'])
    M_dm = sum(simdm['m'])
    SFR = sum(simstar['m'].values[simstar['age'].values < 100]
              ) / 100e6           # Msun/yr
    SFRsd = SFR / (np.pi * R_gal**2.)

    try:
        Zsfr = sum(simgas['Z'] * simgas['SFR']) / sum(simgas['SFR'])
    except ZeroDivisionError:
        # SFR is too small
        Zsfr = np.nan

    # Print properties
    sep = ['30', '10', '20', '40']
    sep1 = '+%' + sep[0] + 's+%' + sep[1] + \
        's+%' + sep[2] + 's+%' + sep[3] + 's+'
    sep2 = '|%' + sep[0] + 's|%' + sep[1] + \
        's|%' + sep[2] + 's|%' + sep[3] + 's|'
    print(sep1 % ((int(sep[0]) * '-'), (int(sep[1]) * '-'),
                  (int(sep[2]) * '-'), (int(sep[3]) * '-')))
    print(sep2 % ('Parameter'.center(int(sep[0])), 'Value'.center(int(
        sep[1])), 'Name in code'.center(int(sep[2])), 'Explanation'.center(int(sep[3]))))
    print(sep1 % ((int(sep[0]) * '-'), (int(sep[1]) * '-'),
                  (int(sep[2]) * '-'), (int(sep[3]) * '-')))
    print(sep2 % ('Stellar mass [1e9 M_sun]'.center(int(sep[0])), str.format(
        "{0:.3f}", M_star / 1e9), "prop['M_star']".center(int(sep[2])), 'Mass of all stellar particles'.center(int(sep[3]))))
    print(sep2 % ('Gas mass [1e9 M_sun]'.center(int(sep[0])), str.format(
        "{0:.3f}", M_gas / 1e9), "prop['M_gas']".center(int(sep[2])), 'Mass of all gas particles'.center(int(sep[3]))))
    print(sep2 % ('SFR [M_sun/yr]'.center(int(sep[0])), str.format("{0:.3f}", SFR), "prop['SFR']".center(
        int(sep[2])), 'SFR averaged over past 100 Myr'.center(int(sep[3]))))
    print(sep2 % ('SFRd [M_sun/yr/kpc^2]'.center(int(sep[0])), str.format(
        "{0:.4f}", SFRsd), "prop['SFRsd']".center(int(sep[2])), 'Surface density of SFR'.center(int(sep[3]))))
    print(sep2 % ('Z [Z_sun]'.center(int(sep[0])), str.format("{0:.3f}", Zsfr), "prop['Z']".center(
        int(sep[2])), 'Mass-weighted metallicity'.center(int(sep[3]))))
    print(sep1 % ((int(sep[0]) * '-'), (int(sep[1]) * '-'),
                  (int(sep[2]) * '-'), (int(sep[3]) * '-')))

    # pdb.set_trace()
    return simgas, simstar, simdm



def simba_to_pd(galnames, raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile, d_data, startGAL=None, endGAL=None, debugplotgas=False, R_max=20.0,
    verbose=False, debug=False):

    """

    Write out useful fields gas, stars, DM particles of the selected galaxies into DataFrames

    Parameters
    ----------
    galnames: str or list of string
        if it's a simple str, it's a file that stores the output from return from select_SFgal_from_simba()
        if it's a list of str, it's output from select_SFgal_from_simba()
    raw_sim_dir: string
        where to look for raw snapshot files
    raw_sim_name_prefix: string
        prefix to the speciifc volume and resolution of snapshots we are looking
    caesar_dir: str
        path to the caesar output files
    name_prefix: str
        predix to file of caesar outputs
    redshiftFile: str
        where to look for the redshift info for this particular set of snapshots
    d_data: str
        directory location to save dataframe files
    startGAL: None or int
        snapshot number to begin extraction. Use to pseudo-parallelized things
        startGAL <= GAL < endGAL
    endGAL: None or int
        snapshot number to end extraction. Use to pseudo-parallelized things
    plotgas: bool
        whether or not to plot gas distribution of galaxy. Uesful for e.g., looking at the galaxy position and see if it might be the COM of two overlapping galaxies - see Romeel's email on Apr 4th, 2019.
    R_max: float
        physical kpc for projection plot
    Returns
    -------
    galnames_selected: list of string
        galaxies names indicating the halo ID, snapshot number, and galaxy ID (numbered based on some predefined criterion)
    zreds_selected: numpy array of float
        redshifts of galaxies


    """

    try:
        import cPickle as pickle
    except:
        import _pickle as pickle
    import caesar
    import yt
    import numpy as np
    import pandas as pd
    import os

    kpc2m = 3.085677580666e19
    resort = False

    # Save the names and redshift for the galaxies that we finally decide to save in DataFrames:
    galnames_selected   =   []
    zreds_selected      =   np.array([])

    if type(galnames) is str:
        try:
            with open(galnames, 'r') as f:
                galnames = pickle.load(f)
        except:
            pass

    try:
        print(galnames['halo'].values, galnames['snap'].values, galnames['GAL'].values, galnames['SFR'].values)
    except TypeError:
        galnames = galnames[0]

    if (sorted(galnames['snap'].values) == galnames['snap'].values).all() == False:
        # then we have to sort galnames by 'snap'
        galnames = galnames.sort_values(['snap', 'SFR'], ascending=[False, False]).reset_index(drop=False)   # drop: avoid old index being added as a column
        resort = True

    snap_hold = []
    for num, (halo, snap, GAL, sfr) in enumerate(zip(galnames['halo'].values, galnames['snap'].values, galnames['GAL'].values, galnames['SFR'].values)):

        if num == 0:
            snap_hold = snap
            s, h, obj, ds, zred = define_ds(snap, raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile)
        else:
            if snap != snap_hold:
                s, h, obj, ds, zred = define_ds(snap, raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile)
                snap_hold = snap

        # if we have sorted galnames by 'snap', we have to do so for obj too.
        if resort:
            gal = obj.galaxies
            obj.galaxies = [gal[i] for i in galnames.index.values]

        if startGAL is None:
            startGAL = 0
        else:
            startGAL = int(startGAL)

        if endGAL is None:
            endGAL = len(galnames.index)
        else:
            endGAL = int(endGAL)

        if startGAL <= GAL < endGAL:
            print('\nNow looking at galaxy # %s with parent halo ID %s in snapshot %s at z = %s' % (
                int(GAL), int(halo), int(snap), zred))

            print('Creating galaxy name:')
            galname = 'h' + str(int(halo)) + '_s' + \
                str(int(snap)) + '_G' + str(int(GAL))
            print(galname)

            galaxy = obj.galaxies[GAL]

            # Exclude galaxies based on global properties to speed things up. TO BE IMPLEMENTED

            if galaxy.masses['gas'] == 0.0 or galaxy.masses['stellar'] == 0.0:
                # then don't bother getting any info on this galaxy
                if galaxy.masses['gas'] == 0.0:
                    print("Don't bother with this galaxy, it has no gas mass.")
                elif galaxy.masses['stellar'] == 0.0:
                    print("Don't bother with this galaxy, it has no stellar mass.")
                continue

            simgas_path = (d_data + 'particle_data/sim_data/z' +
                                 '{:.2f}').format(float(zred)) + '_' + galname + '_sim.gas'
            simstar_path = (d_data + 'particle_data/sim_data/z' +
                                      '{:.2f}').format(float(zred)) + '_' + galname + '_sim.star'
            simdm_path = (d_data + 'particle_data/sim_data/z' +
                                    '{:.2f}').format(float(zred)) + '_' + galname + '_sim.dm'

            if not os.path.exists(simgas_path) and not os.path.exists(simstar_path) and not os.path.exists(simdm_path):

                # Get location and radius for each galaxy belonging to this haloID:
                loc = galaxy.pos            # .in_units('unitary')
                R_gal = galaxy.radius       # kpccm, i.e., co-moving
                # print(galaxy.radii)
                print('Cut out a sphere with radius %s, %s' % (R_gal, R_gal.in_units('kpc')))
                sphere = ds.sphere(loc, R_gal)

                gas_pos = sphere['PartType0', 'Coordinates'].in_units('kpc')
                print('%s SPH particles' % len(gas_pos))

                star_pos_all = sphere['PartType4', 'Coordinates'].in_units('kpc')
                print('%s Star particles' % len(star_pos_all))

                if debug:
                    print("List all stuff inside the raw sim .hdf5")
                    rawSim = raw_sim_dir + raw_sim_name_prefix + \
                            '{:>03}'.format(int(snap)) + '.hdf5'
                    os.system('h5ls -r ' + rawSim)

                    print("")
                    print(ds.field_list)

                if len(gas_pos) > 0 and len(star_pos_all) > 0:
                    print('Extracting all gas particle properties...')
                    gas_pos = gas_pos - loc

                    # rotation to ensure that the galaxies lie in the xy-plane. That makes it easier to visualize later, and you can always add random viewing angles later as well.
                    gas_pos = caesar.utils.rotator(gas_pos,
                                                   galaxy.rotation_angles['ALPHA'],
                                                   galaxy.rotation_angles['BETA'])
                    gas_posx, gas_posy, gas_posz = gas_pos[
                        :, 0].d, gas_pos[:, 1].d, gas_pos[:, 2].d

                    gas_vel = sphere['PartType0', 'Velocities'].in_cgs() / 1.e5
                    gas_vel = caesar.utils.rotator(gas_vel, galaxy.rotation_angles[
                                                   'ALPHA'], galaxy.rotation_angles['BETA'])
                    gas_velx, gas_vely, gas_velz = gas_vel[
                        :, 0].d, gas_vel[:, 1].d, gas_vel[:, 2].d

                    print("gas_velx length: ", len(gas_velx))
                    gas_densities = sphere['PartType0', 'Density'].in_cgs()
                    gas_f_H2 = sphere['PartType0', 'FractionH2']
                    gas_f_neu = sphere['PartType0', 'NeutralHydrogenAbundance']
                    gas_m = sphere['PartType0', 'Masses'].in_units('Msun')
                    # electrons per Hydrogen atom (max: 1.15)
                    gas_x_e = sphere['PartType0', 'ElectronAbundance']

                    # ionized gas mass fraction (because we don't trust
                    # NeutralHydrogenAbundance -- which only includes atomic gas but not moleuclar)
                    # At the end, calculated from
                    gas_f_ion = gas_x_e / max(gas_x_e)
                    gas_f_HI = 1 - gas_f_ion

                    print('\nChecking molecular gas mass fraction from simulation:')
                    print('%.3s %% \n' %
                          (np.sum(gas_m * gas_f_H2) / np.sum(gas_m) * 100.))

                    # Tk
                    gas_Tk = sphere['PartType0', 'Temperature']
                    gas_h = sphere['PartType0', 'SmoothingLength'].in_units(
                        'kpc')       # Tk
                    # Msun/yr
                    gas_SFR = sphere['PartType0', 'StarFormationRate']
                    # gas_SFR = galaxy.sfr
                    gas_Z = sphere['PartType0', 'Metallicity_00'].d / \
                        0.0134               # from RD
                    gas_a_He = sphere['PartType0', 'Metallicity_01'].d
                    gas_a_C = sphere['PartType0', 'Metallicity_02'].d
                    gas_a_N = sphere['PartType0', 'Metallicity_03'].d
                    gas_a_O = sphere['PartType0', 'Metallicity_04'].d
                    gas_a_Ne = sphere['PartType0', 'Metallicity_05'].d
                    gas_a_Mg = sphere['PartType0', 'Metallicity_06'].d
                    gas_a_Si = sphere['PartType0', 'Metallicity_07'].d
                    gas_a_S = sphere['PartType0', 'Metallicity_08'].d
                    gas_a_Ca = sphere['PartType0', 'Metallicity_09'].d
                    gas_a_Fe = sphere['PartType0', 'Metallicity_10'].d

                    print('Extracting all star particle properties...')
                    star_pos = star_pos_all - loc
                    star_pos = caesar.utils.rotator(star_pos, galaxy.rotation_angles[
                                                    'ALPHA'], galaxy.rotation_angles['BETA'])
                    star_posx, star_posy, star_posz = star_pos[
                        :, 0].d, star_pos[:, 1].d, star_pos[:, 2].d

                    star_vel = sphere['PartType4', 'Velocities'].in_cgs() / 1e5
                    star_vel = caesar.utils.rotator(star_vel, galaxy.rotation_angles[
                                                    'ALPHA'], galaxy.rotation_angles['BETA'])
                    star_velx, star_vely, star_velz = star_vel[
                        :, 0].d, star_vel[:, 1].d, star_vel[:, 2].d

                    star_m = sphere['PartType4', 'Masses'].in_units('Msun')
                    # star_m = galaxy.masses['stellar']
                    star_a_C = sphere['PartType4', 'Metallicity_02'].d
                    star_a_O = sphere['PartType4', 'Metallicity_04'].d
                    star_a_Si = sphere['PartType4', 'Metallicity_07'].d
                    star_a_Fe = sphere['PartType4', 'Metallicity_10'].d
                    star_Z = sphere['PartType4', 'Metallicity_00'].d / \
                        0.0134               # from RD

                    print('Info for this snapshot:')
                    omega_baryon = obj.simulation.omega_baryon
                    omega_matter = obj.simulation.omega_matter
                    hubble_constant = obj.simulation.hubble_constant
                    print('Omega baryon: %s' % omega_baryon)
                    print('Omega matter: %s' % omega_matter)
                    print('Hubble constant: %s' % hubble_constant)
                    print('XH: %s' % obj.simulation.XH)

                    current_time = ds.current_time.in_units('yr') / 1.e6  # Myr
                    # in scale factors, do as with Illustris
                    star_formation_a = sphere['PartType4', 'StellarFormationTime'].d
                    star_formation_z = 1. / star_formation_a - 1
                    # Code from yt project (yt.utilities.cosmology)
                    star_formation_t = 2.0 / 3.0 / np.sqrt(1 - omega_matter) * np.arcsinh(np.sqrt(
                        (1 - omega_matter) / omega_matter) / np.power(1 + star_formation_z, 1.5)) / (hubble_constant)  # Mpc*s/(100*km)
                    star_formation_t *=  kpc2m / 100. / (1e6 * 365.25 * 86400)  # Myr
                    star_age = current_time.d - star_formation_t

                    print('Extracting all DM particle properties...')
                    dm_pos_all = sphere['PartType1', 'Coordinates'].in_units('kpc')
                    dm_pos = dm_pos_all - loc
                    dm_pos = caesar.utils.rotator(dm_pos, galaxy.rotation_angles[
                                                  'ALPHA'], galaxy.rotation_angles['BETA'])
                    dm_posx, dm_posy, dm_posz = dm_pos[:, 0].d, dm_pos[:, 1].d, dm_pos[:, 2].d
                    dm_vel = sphere['PartType1', 'Velocities'].in_cgs() / 1e5
                    dm_vel = caesar.utils.rotator(dm_vel, galaxy.rotation_angles[
                                                  'ALPHA'], galaxy.rotation_angles['BETA'])
                    dm_velx, dm_vely, dm_velz = dm_vel[:, 0].d, dm_vel[:, 1].d, dm_vel[:, 2].d
                    dm_m = sphere['PartType1', 'Masses'].in_units('Msun')


                    print('Extracting all BH particle properties...')
                    BH_pos = sphere['PartType5', 'Coordinates'].in_units('kpc') - loc
                    print('%s BH particles' % len(BH_pos))

                    if len(BH_pos) > 0:
                        BH_pos = caesar.utils.rotator(BH_pos,
                                                       galaxy.rotation_angles['ALPHA'],
                                                       galaxy.rotation_angles['BETA'])
                        BH_posx, BH_posy, BH_posz = BH_pos[
                            :, 0].d, BH_pos[:, 1].d, BH_pos[:, 2].d

                        # BH mass
                        # which grows through accretion and mergers w/ other BHs
                        BH_m = sphere['PartType5', 'BH_Mass'].in_units('Msun')

                        # dynamical mass, which enters into gravity calculation
                        BH_m2 = sphere['PartType5', 'Masses'].in_units('Msun')

                        BH_mdot = sphere['PartType5', 'BH_Mdot']

                        # the one that actually matters is the most massive BH particle
                        idx = np.where(BH_m == BH_m.max())
                        BH_m = BH_m[idx]
                        BH_mdot = BH_mdot[idx]

                        print(BH_m, BH_mdot)

                        simbh = pd.DataFrame({'mBH_Msun': BH_m,
                                                'mdot': BH_mdot})

                        simbh.to_pickle(d_data + 'particle_data/sim_data/z' +
                                         '{:.2f}'.format(float(zred)) + '_' + galname + '_sim.bh')

                    # Put into dataframes:
                    simgas = pd.DataFrame({'x': gas_posx, 'y': gas_posy, 'z': gas_posz,
                                           'vx': gas_velx, 'vy': gas_vely, 'vz': gas_velz,
                                           'SFR': gas_SFR, 'Z': gas_Z, 'nH': gas_densities, 'Tk': gas_Tk, 'h': gas_h,
                                           'f_HI1': gas_f_HI,    # atomic and molecular H
                                           'f_neu': gas_f_neu,   # atomic H
                                           'f_H21': gas_f_H2, 'm': gas_m,
                                           'a_He': gas_a_He, 'a_C': gas_a_C, 'a_N': gas_a_N, 'a_O': gas_a_O, 'a_Ne': gas_a_Ne, 'a_Mg': gas_a_Mg,
                                           'a_Si': gas_a_Si, 'a_S': gas_a_S, 'a_Ca': gas_a_Ca, 'a_Fe': gas_a_Fe})
                    simstar = pd.DataFrame({'x': star_posx, 'y': star_posy, 'z': star_posz,
                                            'vx': star_velx, 'vy': star_vely, 'vz': star_velz,
                                            'Z': star_Z, 'm': star_m, 'age': star_age})

                    simdm = pd.DataFrame({'x': dm_posx, 'y': dm_posy, 'z': dm_posz,
                                          'vx': dm_velx, 'vy': dm_vely, 'vz': dm_velz, 'm': dm_m})

                    # Center in position and velocity
                    simgas, simstar, simdm = center_cut_galaxy(
                        simgas, simstar, simdm, plot=False)

                    # Calculate approximate 100 Myr average of SFR for this galaxy
                    m_star = simstar['m'].values
                    age_star = simstar['age'].values
                    SFR_avg = sum(m_star[age_star < 100]) / 100.e6
                    print('Estimated SFR for this galaxy averaged over past 100 Myr: ' + str(SFR_avg))
                    print('SFR in simulation: %s' % galaxy.sfr)
                    print('SFR for parent halo: %s' % galaxy.halo.sfr)

                    print('Stellar mass: %s Msun' % np.sum(m_star))
                    print('Dark matter mass: %s ' % np.sum(dm_m))

                    print('Saving galaxy data as DataFrames in {}').format(d_data)
                    # replace w/ this line if we decide to merge this func to load_module.py and call via parameters_zx.txt
                    # print('Saving galaxy data as DataFrames in {}, as defined in __init__.py').format(d_data)

                    savepath = d_data + 'particle_data/sim_data/'
                    if not os.path.exists(savepath):
                        os.makedirs(savepath)

                    simgas.to_pickle(simgas_path)
                    simstar.to_pickle(simstar_path)
                    simdm.to_pickle(simdm_path)

                    if plotgas:
                        import matplotlib.pyplot as plt
                        from mpl_toolkits.axes_grid1 import AxesGrid
                        # projection plot
                        plt.close('all')

                        savepath = 'plots/sims/'
                        if not os.path.exists(savepath):
                            os.makedirs(savepath)
                        ppp = yt.ProjectionPlot(ds, 0, [('gas', 'density')],
                                             center=sphere.center.value,
                                             width=(R_max, 'kpc'),
                                               # center='c',
                                             weight_field='density')
                        try:
                            ppp.annotate_timestamp(corner='upper_left',
                                                redshift=True,
                                                time=False,
                                                draw_inset_box=True)
                            #p.annotate_scale(corner='upper_right')
                            ppp.annotate_particles((R_max, 'kpc'),
                                                  p_size=20,
                                                  ptype='PartType5',
                                                  minimum_mass=1.e7)
                        except:
                            pass

                        filename = 'z' + '{:.2f}'.format(float(zred)) + '_' + galname + '_gas.pdf'
                        p.save(os.path.join(savepath, filename))

                elif len(gas_pos) == 0 or len(star_pos_all) == 0:

                    if len(star_pos_all) == 0:
                        print("Not saving this galaxy {:} because it has no stellar mass...").format(galname)
                    elif len(gas_pos) == 0:
                        print("Not saving this galaxy {:} because it has no gas mass...").format(galname)

                    if debug:
                        import matplotlib.pyplot as plt
                        from mpl_toolkits.axes_grid1 import AxesGrid
                        # projection plot
                        plt.close('all')

                        savepath = 'plots/sims/'
                        if not os.path.exists(savepath):
                            os.makedirs(savepath)

                        filename = 'z' + '{:.2f}'.format(float(zred)) + '_' + galname + '_gas.pdf'

                        if not os.path.isfile(filename):
                            print("Plotting gas particles...")
                            ppp = yt.ProjectionPlot(ds, 0, [('gas', 'density')],
                                                 center=sphere.center.value,
                                                 width=(R_max, 'kpc'),
                                                   # center='c',
                                                 weight_field='density')
                            try:
                                ppp.annotate_timestamp(corner='upper_left',
                                                    redshift=True,
                                                    time=False,
                                                    draw_inset_box=True)
                                #p.annotate_scale(corner='upper_right')
                                ppp.annotate_particles((R_max, 'kpc'),
                                                      p_size=20,
                                                      ptype='PartType5',
                                                      minimum_mass=1.e7)
                            except:
                                pass

                            ppp.save(os.path.join(savepath, filename))

            else:
                print("Skipping... Already extracted...")

            galnames_selected.append(galname)
            zreds_selected = np.append(zreds_selected, float(zred))

    return galnames_selected, zreds_selected


def pd_bookkeeping(galnames_selected, zreds_selected, zCloudy):
    try:
        import cPickle
    except:
        import _pickle as cPickle
    import os
    models = {'galnames_unsorted': galnames_selected,
              'zreds_unsorted': zreds_selected}
    # call by global_results.py
    outname = '/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/galaxies/z' + str(int(zCloudy)) + '_extracted_gals'
    if os.path.exists(outname):     # make a back up copy if exist
        os.system('mv ' + outname + ' ' + outname + '.bak')
    cPickle.dump(models, open(outname, 'wb'))
    print('Number of galaxies in entire sample extracted from Simba data to : {}'.format(str(len(galnames_selected))))

    return galnames_selected, zreds_selected


def fetch_ss_from_closest_redshift(redshift, redshiftFile):
    """
    Given redshift, fetch the closest snapshot number
    """
    import numpy as np

    _, zs_table, snaps_table = np.loadtxt(redshiftFile, unpack=True)

    snapshotNum = np.argmin(abs(zs_table - redshift))
    return snapshotNum


def hack_galNames(pdPath):
    """

    Since it takes time to extract galaxies, I want to be able to look at properties of galaxies we have extracted so far. So I needed a way to "re-created" the output variables of simba_to_pd().

    This is a very hacky, ad-hoc code. Don't use for science!

    pdPath: str
        path to where the pandas dataframes are saved.
        .gas, .star, .dm
        default should be inside d_data + 'particle_data/sim_data/'

    """

    import glob, os
    gasF = glob.glob(pdPath + "*.gas")
    gasF = map(os.path.basename, gasF)

    zreds_selected = [float(ii[1:5]) for ii in gasF]
    print(zreds_selected)

    galnames_selected = [ii[6:ii.find('_sim')] for ii in gasF]
    print(galnames_selected)

    return galnames_selected, zreds_selected


if __name__ == '__main__':

    raw_sim_dir = '/disk01/rad/sim/m25n256/s48/'
    raw_sim_name_prefix = 'snap_m25n256_'
    caesar_dir = '/disk01/rad/sim/m25n256/s48/Groups/'
    name_prefix = 'm25n256_'
    redshiftFile = '/home/rad/gizmo-extra/outputs_boxspace50.info'

    snapRange = [150, 151]       # snaptable
    zCloudy = 0
    d_data = '/home/dleung/Downloads/SIGAME_dev/sigame/temp/z' + str(int(zCloudy)) + '_data_files/'  #z0

    Nhalo = 2  # 10
    Ngalaxies = 2    # 20

    get_basic_info_from_caesarCat(snapRange, Nhalo, Ngalaxies, caesar_dir, name_prefix)
    ggg = select_SFgal_from_simba(raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, snapRange, Ngalaxies)

    xx, yy = simba_to_pd(ggg, raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile, d_data)

