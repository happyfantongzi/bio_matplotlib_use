#! /usr/bin/env python3
import re,sys,os
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')
from ComputMatrix import subMatrix
import matplotlib.gridspec as gridspec

'''
	this is a temp use script ,the target of this script is to test and find the disadvantage of Computer matrix module
	and update 'The focus region'
	so ,how to do it?
		1.collect matrix files
		2.use Gridspec to layout
		3.focus on detail like colorbar,the fontsize
		4.plot the whole genome then part of it
		5.should I add other information in it? may be i will add boundary information
'''
#多细胞系路径:  /annoroad/data1/bioinfo/PMO/Public/database/HIC/BJ161163-03/data/Fire_data/GSE87112/Pre_deal/contact_maps/contact_maps/ICE/process/process/EdTAD/
#chr19 0 61342430  resolution=40000 
list_cells_df=pd.read_table("/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/process/multi_cells_TAD/list.cells",names=['name','path'])
chrom_size_file="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/process/multi_cells_TAD/chrom_mm9.sizes"
chrom_size_df=pd.read_table(chrom_size_file,names=['chr','size'])
for ii ,rr in chrom_size_df.iterrows():
	out_fig="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/process/multi_cells_TAD2/{}_cell_tad.pdf".format(rr['chr'])
	print(out_fig)
	fig=plt.figure(figsize=(16,4))
	G=gridspec.GridSpec(list_cells_df.shape[0],20)
	plt.suptitle('{} TAD'.format(rr['chr']),fontsize=15)
	for i,r in list_cells_df.iterrows():
		mat=subMatrix(rr['chr'],0,rr['size'],40000,0)
		mat.sample=r['name']
		mat_path=r['path'].replace('chr19',rr['chr'])
		print(rr['chr'])
		print(mat_path)
		mat.mat_path=mat_path
		mat.colormap=plt.get_cmap('Reds')
		ax = plt.subplot(G[i:i+1,:],axisbg='white')
		ax.set_ylabel(r['name'],rotation='horizontal',horizontalalignment='right')
		mat.tri_map(ax,mat.all_mat,False)
	fig.savefig(out_fig)

