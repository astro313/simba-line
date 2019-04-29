"""

getting all z~6 galaxies in the 1024 sims, from the 25 cMpc/h box, but they seem to have too low SFR - so for my purpose, I maybe should look into the 50 box instead.

"""

import sys
if sys.version_info[0] < 3:
    raise Exception("Python 3 or a more recent version is required.")

import _pickle as pickle     # py3
from parse_simba import *
import socket
import os

if len(sys.argv) > 1:
    debug = sys.argv[1]
else:
    debug = False

snapRange = [36]    # don't put 036!
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

# Nhalo and Ngalaxies below don't affect the ouput, just for inspection
# total_Ngal = get_basic_info_from_caesarCat(snapRange, Nhalo=10, Ngalaxies=10, caesar_dir, name_prefix)
# print "Using all galaxies found from the Caesar catalog"
# print total_Ngal

# ggg = select_SFgal_from_simba(raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, snapRange, Ngalaxies=total_Ngal, saveggg='ggg.pkl')

galnames, zred = simba_to_pd('ggg.pkl', raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile, d_data, debug=debug, startGAL=None, endGAL=None)

_, _ = pd_bookkeeping(galnames, zred, zCloudy)


# to go with whatever has been extracted so far.... then uncomment the following...
# from parse_simba import hack_galNames
# galnames, zred = hack_galNames(d_data+'particle_data/sim_data/')
# _, _ = pd_bookkeeping(galnames, zred, zCloudy)

