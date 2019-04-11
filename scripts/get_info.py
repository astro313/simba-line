import caesar
import numpy as np
import cPickle as pickle
import socket
import os
import sys
from readgadget import readsnap

def info(obj, snapFile, top=None, savetxt=False):
    """
    Customized version to print more info on the CAESAR extracted galaxies.

    Parameters
    ----------
    obj : :class:`main.CAESAR`
        Main CAESAR object.
    top : int
        Number of objects to print.

    """

    group_list = obj.galaxies
    nobjs = len(group_list)

    if top > nobjs or top is None:
        top = nobjs

    if obj.simulation.cosmological_simulation:
        time = 'z=%0.3f' % obj.simulation.redshift
    else:
        time = 't=%0.3f' % obj.simulation.time

    output  = '\n'
    output += '## Largest %d Galaxies \n' % (top)
    if hasattr(obj, 'data_file'): output += '## from: %s\n' % obj.data_file
    output += '## %d @ %s' % (nobjs, time)
    output += '\n\n'

    cnt = 1

    output += ' ID      Mstar     Mgas      MBH    fedd    SFR [Msun/yr]       r_baryon   r_gas      r_gas_half_mass      r_stellar    r_stellar_half_mass    Z_sfrWeighted [/Zsun]    Z_massWeighted     Z_stellar     T_gas_massWeighted    T_gas_SFRWeighted   fgas   nrho      Central\t|  Mhalo     HID\n'
    output += ' ----------------------------------------------------------------------------------------\n'

    h = obj.simulation.hubble_constant
    bhmdot = readsnap(snapFile, 'BH_Mdot','bndry',suppress=1)*1.e10/h/3.08568e+16*3.155e7 # in Mo/yr

    for ii, o in enumerate(group_list):
        phm, phid = -1, -1

        if o.halo is not None:
            phm, phid = o.halo.masses['total'], o.halo.GroupID

        try:
            bhmdots = [bhmdot[k] for k in o.bhlist]
            imax = np.argmax(o.masses['bh'])
            bm = o.masses['bh'][imax]

            frad = 0.1
            mdot_edd = 4*np.pi*6.67e-8*1.673e-24/(frad*3.e10*6.65245e-25) * bm * 3.155e7 # in Mo/yr
            bmdot = bhmdots[imax]        # only the massive BH particle matters.
            fedd = bmdot / mdot_edd
        except:
            bm = 0
            fedd = 0

        output += ' %04d  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e   %0.3f   %0.3f  %0.3f  %0.2e  %0.2e  %0.3f  %0.2e  %s\t|  %0.2e  %d \n' % \
                  (o.GroupID, o.masses['stellar'], o.masses['gas'],
                   bm,
                   fedd,
                   o.sfr,
                   o.radii['baryon'],
                   o.radii['gas'],
                   o.radii['gas_half_mass'],
                   o.radii['stellar'],
                   o.radii['stellar_half_mass'],
                   o.metallicities['sfr_weighted']/0.0134,
                   o.metallicities['mass_weighted'],
                   o.metallicities['stellar'],
                   o.temperatures['mass_weighted'],
                   o.temperatures['sfr_weighted'],
                   o.gas_fraction,
                   o.local_number_density, o.central,
                   phm, phid)

        cnt += 1
        if cnt > top: break

    print(output)

    if savetxt:
        outputFileName = snapFile[50:-5] + '_top' + str(top) + '.txt'
        f = open(outputFileName, 'w')
        f.write(output)
        f.close()

    return output


def setup_cmap(cm='magma'):    # "plasma_r"
    import matplotlib.pyplot as plt
    return plt.set_cmap(cm)


def setup_plot():
    import matplotlib
    from cycler import cycler
    # print(matplotlib.matplotlib_fname())

    matplotlib.rcParams.update({'figure.figsize': (8, 5)    # inches
                                , 'font.size': 16      # points
                                , 'legend.fontsize': 10      # points
                                , 'lines.linewidth': 2       # points
                                , 'axes.linewidth': 1       # points
                                , 'axes.prop_cycle': cycler('color', 'bgrcmyk')
                                , 'text.usetex': True
                                , 'font.family': "serif"  # Use serifed fonts
                                , 'xtick.major.size': 13     # length, points
                                , 'xtick.major.width': 1     # points
                                , 'xtick.minor.size': 8     # length, points
                                , 'xtick.minor.width': 1     # points
                                , 'ytick.major.size': 13     # length, points
                                , 'ytick.major.width': 1     # points
                                , 'ytick.minor.size': 8     # length, points
                                , 'xtick.labelsize': 14
                                , 'ytick.labelsize': 14
                                , 'ytick.minor.width': 1     # points
                                , 'font.serif': ("times", "Computer Modern Roman", "New Century Schoolbook", "Bookman"), 'font.sans-serif': ("Helvetica", "Avant Garde", "Computer Modern Sans serif"), 'font.monospace': ("Courier", "Computer Modern Typewriter"), 'font.cursive': "Zapf Chancery"
                                })
    cm = setup_cmap()
    return cm


def plot_info(colNumx, colNumy, inFile, logx=True, logy=True, xlabel='',
              ylabel='',
              tag='',
              saveFig=True):
    """

    colNumx: int
        column to read from the .txt file, used as x-axis in plotting

    colNumy: int
        column to read from the .txt file, used as x-axis in plotting

    inFile: str
        name of input file to read from where we extract global properties of galaxies


    """

    import matplotlib.pyplot as plt

    cm = setup_plot()

    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(111)




if __name__ == '__main__':

    if len(sys.argv) > 1:
        debug = sys.argv[1]
        LoadHalo = False
    else:
        debug = False
        LoadHalo = True

    snapRange = [36]    # don't put 036
    zCloudy = 6
    raw_sim_name_prefix = 'snap_m25n1024_'
    name_prefix = 'm25n1024_'

    host = socket.gethostname()
    if 'ursa' in host:
        raw_sim_dir = '/disk01/rad/sim/m25n1024/s50/'
        caesar_dir = '/disk01/rad/sim/m25n1024/s50/Groups/'
        redshiftFile = '/home/rad/gizmo-extra/outputs_boxspace50.info'
        d_data = '/home/dleung/Downloads/SIGAME_dev/sigame/temp/z' + str(int(zCloudy)) + '_data_files/'
    elif 'flatironinstitute.org' or 'worker' in host:
        raw_sim_dir = '/mnt/ceph/users/daisyleung/simba/sim/m25n1024/s50/'  # dummy
        caesar_dir = '/mnt/ceph/users/daisyleung/simba/sim/m25n1024/s50/Groups/'
        redshiftFile = '/mnt/ceph/users/daisyleung/simba/gizmo-extra/outputs_boxspace50.info'
        d_data = '/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/z' + str(int(zCloudy)) + '_data_files/'

    infile = caesar_dir + name_prefix + '{:0>3}'.format(int(36)) + \
                '.hdf5'
    obj = caesar.load(infile, LoadHalo=LoadHalo)
    snapFile = raw_sim_dir + raw_sim_name_prefix + '{:0>3}'.format(int(36)) + '.hdf5'

    output = info(obj, snapFile, top=None, savetxt=True)

