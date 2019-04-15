"""

getting all z~6 galaxies in the 1024 sims, from the 25 cMpc/h box, but they seem to have too low SFR - so for my purpose, I maybe should look into the 50 box instead.

"""

import cPickle as pickle
try:
    from parse import *
except:
    from parse_simba import *
import socket
import os
import sys

if len(sys.argv) > 1:
    debug = sys.argv[1]
else:
    debug = False

print("Extracting galaxies from z6 ONLY, for now.. we may need to extract galaxies from other redshift.... But hopefully we don't have to")
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

# this doesn't affect ouput, just for inspection
# Nhalo = 10
# Ngalaxies = 10
# total_Ngal = get_basic_info_from_caesarCat(snapRange, Nhalo, Ngalaxies, caesar_dir, name_prefix)

# print "Using all galaxies found from the Caesar catalog"
# print total_Ngal
# Ngalaxies = total_Ngal

# ggg = select_SFgal_from_simba(raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, snapRange, Ngalaxies, saveggg='ggg.pkl')

galnames, zred = simba_to_pd('ggg.pkl', raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile, d_data, plotgas=False, debug=debug)

_, _ = pd_bookkeeping(galnames, zred, zCloudy)


# to go with whatever has been extracted so far.... then uncomment the following...
# from parse_simba import hack_galNames
# galnames, zred = hack_galNames(d_data+'particle_data/sim_data/')
# _, _ = pd_bookkeeping(galnames, zred, zCloudy)

