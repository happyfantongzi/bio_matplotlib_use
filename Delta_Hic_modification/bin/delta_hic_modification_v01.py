#! /usr/bin/env python3
'''
	this module is used to integret the class of 
		1.delta matrix 
		2.histone modification
		3.delta tad boundaries
	and all the focus is below:
		1.xlim of every subplot
		2.the Gridspect of every subplot
		3.draw genes and annotation.....
		4.lay out....
'''
import os,sys,re
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
plt.style.use('ggplot')
Bin = os.path.abspath(os.path.dirname(__file__))
import unittest
from delta_matrix import delta_hic
from multi_modification import multi_modification
from delta_tad_bounds import delta_boundaries

__author__='fanxuning'
__mail__='xuningfan@genome.cn'

class delta_hic_modification:
	def __init__(self,chr,start,end):
		self.chr=chr
		self.start=start
		self.end=end
		self.sam1=''
		self.sam2=''
		self.mat1_path=''
		self.mat2_path=''
		self.bound1_path=''
		self.bound2_path=''
		self.listbams=''
		self.exists_fs_file=None
		self.delta_matrix_grid=5
		self.delta_bound_grids=1
		self.genome_start=0
		self.genome_end=0
		self.delta_Gs=gridspec.GridSpec(10,1)
		self.layouts=[]
		self.outdir=''
	
	@property
	def multi_modification_grids(self):
		if len(self.exists_fs_file)==0:
			return len(self.exists_fs_file)
		else:
			listbam_df=pd.read_table(self.listbams,names=['names','treat','control'])
			return (listbam_df.shape[0])

	def delta_matrix_obj(self,ax):
		dh=delta_hic(self.chr,self.start,self.end)
		dh.sam1=self.sam1
		dh.sam2=self.sam2
		dh.mat1_path=self.mat1_path
		dh.mat2_path=self.mat2_path
		raw_mat=dh.minus_matrix()
		dh.draw_minus_matrix(ax,raw_mat)
		return dh
	
	def delta_bds_obj(self,ax):  
		debs=delta_boundaries()
		debs.bound1_path=self.bound1_path
		debs.bound2_path=self.bound2_path
		debs.chr=self.chr
		debs.start=self.start
		debs.end=self.end
		debs.classfied_boundaries()
		debs.draw_classify_domain(ax)
		self.genome_start=debs.genome_start
		self.genome_end=debs.genome_end
		return debs

	def multi_modfication_obj(self,layouts):
		modi=multi_modification()
		modi.chr=self.chr
		print('modification start end...')
		modi.mstart=self.genome_start
		modi.mend=self.genome_end
		print('modification start end:{}\t{}'.format(modi.mstart,modi.mend))
		modi.bams_list=self.listbams
		modi.outdir=self.outdir
		#modi.multi_processing()
		modi.exists_fs_file=self.exists_fs_file
		modi.gs=self.delta_Gs
		modi.modification_layout=layouts
		return modi

	def pre_delta_lay_out(self):
		'''
			Three  subunits 
				1.delta Hic matrix :   4 grids  layouts=[0:3,:]
				2.delta tad bound  :   2 grids  layouts=[3:4,:]
				3.multi modification : 4 grids  layouts=[4:,:]
		'''
		top_grids=self.delta_matrix_grid+self.delta_bound_grids+2
		grids=top_grids+self.multi_modification_grids
		self.delta_Gs=gridspec.GridSpec(grids,1)
		layouts=list(range(top_grids,grids))
		return layouts

	def draw_all_delta_map(self,fig):
		Chip_layouts=self.pre_delta_lay_out()
		G=self.delta_Gs
		#--------------draw delta matrix -------------------------#
		gridsSub1=gridspec.GridSpecFromSubplotSpec(1,1,
				subplot_spec=G[0:self.delta_matrix_grid,:],wspace=0.5,hspace=0.7)
		ax1=fig.add_subplot(gridsSub1[:,:],axisbg='white')
		self.delta_matrix_obj(ax1)
		#---------------draw delta tad ---------------------------#
		gridsSub2=gridspec.GridSpecFromSubplotSpec(1,1,
				subplot_spec=G[self.delta_matrix_grid+2:(self.delta_matrix_grid+self.delta_bound_grids+2),:],wspace=0.5,hspace=0.7)
		ax2=fig.add_subplot(gridsSub2[:,:],axisbg='white')
		self.delta_bds_obj(ax2)
		#--------------draw multi_modfication---------------------#
		self.multi_modfication_obj(Chip_layouts).multi_draw_ChIP_signal(fig)
		print ("all map....")

class test_delta_hic_modification(unittest.TestCase):
	def test_draw(self):
		dhm=delta_hic_modification('chr19',10000000,48342430)
		data_dir="{}/../Example/Delta_Hic_modification/data/".format(Bin)
		control_path='{}/GM.mm9_chr19.mm9_chr19.matrix'.format(data_dir)
		treat_path='{}/DM.mm9_chr19.mm9_chr19.matrix'.format(data_dir)
		dhm.sam1='GM';dhm.mat1_path=control_path
		dhm.sam2='DM';dhm.mat2_path=treat_path
		dhm.bound1_path="{}/GM_chr19.boundaries".format(data_dir)
		dhm.bound2_path="{}/DM_chr19.boundaries".format(data_dir)
		listbams="{}/bamfile.list".format(data_dir)
		dhm.listbams=listbams
		hists_fs_dir="{}/CHIP_Signal/".format(data_dir)
		hists=['GH3K27ac_class_data.fs','DH3K27ac_class_data.fs','GH3K27me3_class_data.fs','DH3K27me3_class_data.fs','GH3K4me1_class_data.fs','DH3K4me1_class_data.fs','GH3K9me3_class_data.fs','DH3K9me3_class_data.fs']
		print('fs_files.....')
		fs_files=[hists_fs_dir+x for x in hists]
		print(fs_files)
		dhm.exists_fs_file=fs_files
		outdir="{}/../Example/Delta_Hic_modification/".format(Bin)
		data_dir="{}/data/CHIP_Signal/".format(outdir)
		dhm.outdir=data_dir
		outfig="{}/delta_hic_modification.pdf".format(outdir)
		print(outfig)
		fig=plt.figure(figsize=(15,2))
		dhm.draw_all_delta_map(fig)
		fig.savefig(outfig)

if __name__=='__main__':
	unittest.main()
