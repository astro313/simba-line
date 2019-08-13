from parse_simba import pd_bookkeeping
import os
import pandas as pd
from yt2caesar import *

zCloudy = 6

if not os.path.exists('xxx50/'):
    os.mkdir('xxx50/')


pp = particles2pd(snapRange=[36],name_prefix='m50n512_', feedback='s50/', zCloudy=zCloudy, user='Daisy', part_threshold=int(64/8), sfr_threshold=0.1, denseGasThres=1.e4)
ggg2, zred = pp.run(savepath='xxx50/', outname='/mnt/home/daisyleung/Downloads/SIGAME_dev/sigame/temp/galaxies/z' + str(int(zCloudy)) + '_extracted_gals_m50_n512', emptyDM=True, caesarRotate=False, LoadHalo=True)

print(" **** \n \
1/ Manually copy the pandas DF to sigame sim_data/ \n \
2/ Manually create sigame ISM_data/ \n \
3/ ln paramters_zx.txt file.. \n \
4/ update galaxies/z6_extracted_gals_xxx file path... \n \
5/ make sure not going to over-write global_results/zx_..._abun_abun file \n ***")
