#! /usr/bin/env python3
import re,sys,os
import argparse
import numpy as np
import pandas as pd
from delta_hic_modification import delta_hic_modification
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')

'''
Data collection:
a:data view region:
	dispearing borders
	reforming borders
b: matrix bounds......
		sam1	sam2
	mat	mat1	mat2
	bound	bound1	bound2
	listbam	listbam1	listbam2
	gene	anno	anno
	'''
#step 1 data view region collections
diff_bound_path="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/Analysis_script/TAD_switch_stat/Example/Analysis20180120/Data_view/"
all_dispearing_path="{}/all_disapearing_borders.bed".format(diff_bound_path)
adf=pd.read_table(all_dispearing_path,names=['chr','start','end'])
adf['type']='dispearing'
all_reforming_path="{}/all_reforming_borders.bed".format(diff_bound_path)
arf=pd.read_table(all_reforming_path,names=['chr','start','end'])
arf['type']='reforming'
diff_bound_df=pd.concat([adf,arf])
diff_bound_df=diff_bound_df.reset_index(drop=True)

annote_file="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/bin/Delta_Hic_modification/Example/Delta_Hic_modification/data/DM_GM.anno.xls"

#step 2 hunting matrix path ,bound_path .etc
from glob import glob
samples=['GM','DM']
mats_path="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/Analysis_standard/Analysis/process/TAD_dekker/"
bounds_path="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/Analysis_standard/Analysis/process/TAD_dekker/"
Chip_path="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/bin/Delta_Hic_modification/Example/Delta_Hic_modification/data/"

outdir='/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/process/'

up_extend=5000000;down_extend=5000000
for i,r  in  diff_bound_df.iterrows():
	print('\n\n\n')
	print("#####################################{}\t {}:{}-{} ......####################################################".format(i,r['chr'],r['start'],r['end']))
	start =r['start']-up_extend if r['start']>up_extend else 0
	end=r['end']+ down_extend 
	#pre data file prepared....
	mat1=glob("{0}/{1}/40000/{2}/{1}.mm9_{2}.mm9_{2}.matrix".format(mats_path,samples[0],r['chr']))[0]
	mat2=glob("{0}/{1}/40000/{2}/{1}.mm9_{2}.mm9_{2}.matrix".format(mats_path,samples[1],r['chr']))[0]
	bound1=glob("{0}/{1}/40000/{2}/{1}_{2}--is1000001--nt0.3--ids400001--ss1--immean.insulation.boundaries".format(bounds_path,samples[0],r['chr']))[0]
	bound2=glob("{0}/{1}/40000/{2}/{1}_{2}--is1000001--nt0.3--ids400001--ss1--immean.insulation.boundaries".format(bounds_path,samples[1],r['chr']))[0]
	print('listbams....')
	listbams="{}/bamfile.list".format(Chip_path)
	print(listbams)
	hists_fs_dir="{}/CHIP_Signal_5k/".format(Chip_path)
	print('hists_fs_dir.......')
	print(hists_fs_dir)
	hists=['GH3K27ac_class_data.fs','DH3K27ac_class_data.fs','GH3K27me3_class_data.fs','DH3K27me3_class_data.fs','GH3K4me1_class_data.fs','DH3K4me1_class_data.fs','GH3K9me3_class_data.fs','DH3K9me3_class_data.fs']
	fs_files=[hists_fs_dir+x for x in hists ]
	exists_minus_mats=glob("{}/dispearing/{}_quantile_norm.mat".format(outdir,r['chr']))
	#step3  build delta hic modification instance....
	dhm=delta_hic_modification(r['chr'],start,end)
	dhm.listbams=listbams
	dhm.sam1=samples[0];dhm.sam2=samples[1]
	dhm.mat1_path=mat1;dhm.mat2_path=mat2
	dhm.bound1_path=bound1;dhm.bound2_path=bound2
	if not len(exists_minus_mats)==0:
		dhm.exists_minus_mat=exists_minus_mats[0]
	dhm.Chip_binSize=5000
	flag=True
	for x in fs_files:
		if len(glob(x))==0:
			flag=False
	if flag:
		print('exists_fs_file....')
		print(fs_files)
		dhm.exists_fs_file=fs_files
	dhm.annote_file=annote_file
	dhm.outdir="{}/{}/".format(outdir,r['type'])
	outfig="{}/{}_{}_{}_delta_hic.pdf".format(dhm.outdir,r['chr'],r['start'],r['end'])
	fig=plt.figure(figsize=(4,2))
	dhm.draw_all_delta_map(fig)
	fig.savefig(outfig)
