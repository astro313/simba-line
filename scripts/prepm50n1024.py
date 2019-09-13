from parse_simba import pd_bookkeeping
import os
import pandas as pd
from yt2caesar import *

zCloudy = 6

import sys
ss = sys.argv[1]

outDIR = 'xxx50_n1024_ss' + ss + '/'
if not os.path.exists(outDIR):
    os.mkdir(outDIR)

pp = particles2pd(snapRange=[int(ss)],name_prefix='m50n1024_', feedback='s50/', zCloudy=zCloudy, user='Daisy', part_threshold=64, sfr_threshold=0.1, denseGasThres=1.e4)
ggg2, zred = pp.run(savepath=outDIR, outname='/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/galaxies/z' + str(int(zCloudy)) + '_extracted_gals_m50_n1024_ss' + ss, emptyDM=True, caesarRotate=True, LoadHalo=True)

print(" **** \n \
1/ Manually copy the pandas DF to sigame sim_data/ \n \
2/ Manually create sigame ISM_data/ \n \
3/ ln paramters_zx.txt file.. \n \
4/ update galaxies/z6_extracted_gals_xxx file path... \n \
5/ make sure not going to over-write global_results/zx_..._abun_abun file \n ***")
