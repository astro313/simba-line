"""

getting all z~6 galaxies in the 1024 sims, from the 25 cMpc/h box


"""

import cPickle as pickle
try:
    from parse import *
except:
    from parse_simba import *

raw_sim_dir = '/disk01/rad/sim/m25n1024/s50/'
raw_sim_name_prefix = 'snap_m25n1024_'
caesar_dir = '/disk01/rad/sim/m25n1024/s50/Groups/'
name_prefix = 'm25n1024_'
redshiftFile = '/home/rad/gizmo-extra/outputs_boxspace50.info'

snapRange = [36, 150]    # don't put 036
zCloudy = 6

dd_data = '/home/dleung/Downloads/SIGAME_dev/sigame/temp/z' + str(int(zCloudy)) + '_data_files/'

multiprocessing = True

# this doesn't affect ouput, just for inspection
# Nhalo = 10
# Ngalaxies = 10
# total_Ngal = get_basic_info_from_caesarCat(snapRange, Nhalo, Ngalaxies, caesar_dir, name_prefix)

# print "Using all galaxies found from the Caesar catalog"
# print total_Ngal
# Ngalaxies = total_Ngal

# ggg = select_SFgal_from_simba(raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, snapRange, Ngalaxies)

if multiprocessing:
    # save output to .txt for multiprocessing
    prefix = raw_sim_dir[raw_sim_dir.find('sim/')+4:-1].replace('/', '_')
    save_SFgal_txt(ggg, prefix)
    split_by_snap()

    # multiprocessing
    for niter, iii in enumerate(snapRange):
        if niter == 0:
            # multiprocessing
            gg, zz = simba_to_pd_multiprocessing(4,  raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile, dd_data, snap=iii)
        else:
            for cc in gg:
                gg.append(cc)
            for cc in zz:
                zz.append(cc)
    galnames, redshifts = pd_bookkeeping(gg, zz, zCloudy)

else:
    # # save output ggg to pickle file for debugging purpose,
    # I don't want to have to run select_SFgal_from_simba() again and again
    # with open('ggg.pkl', 'w') s f:
    #     pickle.dump([ggg], f)

    with open('ggg.pkl', 'r') as f:
        ggg = pickle.load(f)

    galnames, zred = simba_to_pd(ggg, raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile, dd_data)




fetch_BH(ggg, raw_sim_dir, raw_sim_name_prefix, caesar_dir, name_prefix, redshiftFile)
