import caesar
import numpy as np
import cPickle as pickle
import socket
import os
import sys
from readgadget import readsnap

def info(obj, snapFile, top, savetxt=False):
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

    if top > nobjs:
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

    output += ' ID      Mstar     Mgas      MBH    MBH_dot [Msun/yr]    SFR [Msun/yr]       r_baryon   r_gas      r_gas_half_mass      r_stellar    r_stellar_half_mass    Z_sfrWeighted [/Zsun]    Z_massWeighted     Z_stellar     T_gas_massWeighted    T_gas_SFRWeighted   fgas   nrho      Central\t|  Mhalo     HID\n'
    output += ' ----------------------------------------------------------------------------------------\n'

    h = obj.simulation.hubble_constant
    bhmdot = readsnap(snapFile, 'BH_Mdot','bndry',suppress=1)*1.e10/h/3.08568e+16*3.155e7 # in Mo/yr

    for ii, o in enumerate(group_list):
        phm, phid = -1, -1

        if o.halo is not None:
            phm, phid = o.halo.masses['total'], o.halo.GroupID

        print o.masses['bh']
        import pdb; pdb.set_trace()

        bhmdots = [bhmdot[k] for k in o.bhlist]
        imax = np.argmax(o.masses['bh'])
        bmdot = bhmdots[imax]        # only the massive BH particle matters.

        print o.masses['bh'], bmdot
        import pdb; pdb.set_trace()

        output += ' %04d  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e  %0.2e   %0.3f   %0.3f  %0.3f  %0.2e  %0.2e  %0.3f  %0.2e  %s\t|  %0.2e  %d \n' % \
                  (o.GroupID, o.masses['stellar'], o.masses['gas'],
                   o.masses['bh'],
                   bmdot,
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
    # if savetxt:
    return output



if __name__ == '__main__':

    if len(sys.argv) > 1:
        debug = sys.argv[1]
    else:
        debug = False

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
    obj = caesar.load(infile)
    snapFile = raw_sim_dir + raw_sim_name_prefix + '{:0>3}'.format(int(36)) + '.hdf5'

    output = info(obj, snapFile, top=1)

