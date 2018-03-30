#! /usr/bin/env python3
import re,sys,os
import pandas as pd
import numpy as np
from delta_matrix import delta_hic
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')

outdir='/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/process/delta_matrix_distr/Exception/'
chrom_size_path="/annoroad/share/software/package/package/hic/HiC-Pro-master/annotation/mm9/genome/process/Prepare/chrom_mm9.sizes"
# /annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/Analysis_standard/Analysis/process/TAD_dekker/DM/40000/chr1/DM.mm9_chr1.mm9_chr1.matrix
chrom_size_df=pd.read_table(chrom_size_path,names=['chr','size'])
for i,r in chrom_size_df.iterrows():
	if r['chr']=='chr1':continue
	print("chr\t{}".format(r['chr']))
	control_path="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/Analysis_standard/Analysis/process/TAD_dekker/GM/40000/{0}/GM.mm9_{0}.mm9_{0}.matrix".format(r['chr'])
	treat_path="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/Analysis_standard/Analysis/process/TAD_dekker/DM/40000/{0}/DM.mm9_{0}.mm9_{0}.matrix".format(r['chr'])
	hic_minus=delta_hic(r['chr'],0,r['size'])
	hic_minus.sam1='GM'
	hic_minus.sam2='DM'
	hic_minus.mat1_path=control_path
	hic_minus.mat2_path=treat_path
	hic_minus.outdir=outdir
	out_boxplot='{}/{}_minus_interaction_dist_boxplot.pdf'.format(outdir,r['chr'])
	print(out_boxplot)
	fig=plt.figure(figsize=(4,6))
	ax=fig.add_axes([0.2,0.1,0.7,0.8],axisbg='white')
	raw_mat=hic_minus.minus_matrix()
	hic_minus.draw_Interaction_Change_distance_distr(ax,0,raw_mat)
	fig.savefig(out_boxplot)

