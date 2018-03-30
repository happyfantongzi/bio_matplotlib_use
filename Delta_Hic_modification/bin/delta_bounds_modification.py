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
from gene_expression import Gene_Express

__author__='fanxuning'
__mail__='xuningfan@genome.cn'

class delta_bounds_modification:
	def __init__(self,chr,start,end):
		self.chr=chr
		self.start=start
		self.end=end
		self.sam1=''
		self.sam2=''
		self.insulat1=''
		self.insulat2=''
		self.bound1_path=''
		self.bound2_path=''
		self.listbams=''
		self.exists_fs_file=None
		self.delta_bound_grids=1
		self.gene_express_grid=2
		self.genome_start=0
		self.genome_end=0
		self.annote_file=''
		self.delta_Gs=gridspec.GridSpec(10,1)
		self.Chip_binSize=10000
		self.Chip_layouts=[]
		self.bounds_layout=[]
		self.outdir=''
	
	@property
	def multi_modification_grids(self):
		if not self.exists_fs_file==None:
			return len(self.exists_fs_file)
		else:
			listbam_df=pd.read_table(self.listbams,names=['names','treat','control'])
			return (listbam_df.shape[0])

	def delta_bds_obj(self,G):  
		debs=delta_boundaries()
		debs.bound1_path=self.bound1_path
		debs.insulat1=self.insulat1
		debs.insulat2=self.insulat2
		debs.bound2_path=self.bound2_path
		debs.sam1=self.sam1
		debs.sam2=self.sam2
		debs.chr=self.chr
		debs.start=self.start
		debs.end=self.end
		debs.classfied_boundaries()
		#debs.draw_classify_domain(ax)
		self.genome_start=debs.genome_start
		self.genome_end=debs.genome_end
		self.delta_bound_grids=debs.pre_delta_bound_layout()
		debs.bi_Gs=G
		return debs

	def multi_modfication_obj(self,layouts):
		modi=multi_modification()
		modi.chr=self.chr
		print('modification start end...')
		modi.mstart=self.genome_start
		modi.mend=self.genome_end
		modi.binSize=self.Chip_binSize
		print('modification start end:{}\t{}'.format(modi.mstart,modi.mend))
		modi.bams_list=self.listbams
		modi.outdir=self.outdir
		modi.exists_fs_file=self.exists_fs_file
		modi.gs=self.delta_Gs
		modi.modification_layout=layouts
		return modi

	@property
	def gene_express(self):
		gene_exp=Gene_Express(self.chr,self.genome_start,self.genome_end)
		gene_exp.annote_file=self.annote_file
		return gene_exp

	def pre_delta_lay_out(self):
		'''
			Three  subunits 
				1.delta tad bound  :   4 grids  layouts=[3:4,:]
				2.multi modification : 4 grids  layouts=[4:,:]
				3.gene_express    : 1 grid
		'''
		self.delta_bound_grids=self.delta_bds_obj(self.delta_Gs).pre_delta_bound_layout()
		top_grids=self.delta_bound_grids+1
		grids=top_grids+self.multi_modification_grids*2+self.gene_express_grid*2+1
		self.delta_Gs=gridspec.GridSpec(grids,1)
		self.Chip_layouts=list(range(top_grids,grids,2))

	def draw_all_delta_map(self,fig):
		self.pre_delta_lay_out()
		Chip_layouts=self.Chip_layouts
		gene_layouts=Chip_layouts[-1]+1
		G=self.delta_Gs
		#---------------draw delta tad ---------------------------#
		self.delta_bds_obj(G).draw_delta_map(fig)
		#--------------draw multi_modfication---------------------#
		self.multi_modfication_obj(Chip_layouts).gs=G
		self.multi_modfication_obj(Chip_layouts).multi_draw_ChIP_signal(fig)

		#--------------draw gene expression-----------------------#
		gridsSub4=gridspec.GridSpecFromSubplotSpec(1,1,
				subplot_spec=G[Chip_layouts[-2]:(Chip_layouts[-1]+self.gene_express_grid),:],wspace=0.5,hspace=0.3)
		ax4=fig.add_subplot(gridsSub4[:,:],axisbg='white')
		self.gene_express.draw_gene_expression(ax4)

		print ("all map....")

class test_delta_hic_modification(unittest.TestCase):
	def test_draw(self):
		dhm=delta_bounds_modification('chr19',10000000,48342430)
		data_dir="{}/../Example/delta_bounds_modification/data/".format(Bin)
		dhm.bound1_path="{}/GM_chr19.boundaries".format(data_dir)
		dhm.bound2_path="{}/DM_chr19.boundaries".format(data_dir)
		dhm.insulat1="{}/GM_chr19.insulation".format(data_dir)
		dhm.insulat2="{}/DM_chr19.insulation".format(data_dir)
		dhm.sam1='GM'
		dhm.sam2='DM'
		dhm.annote_file="{}/DM_GM.anno.xls".format(data_dir)
		listbams="{}/bamfile.list".format(data_dir)
		dhm.listbams=listbams
		hists_fs_dir="{}/CHIP_Signal_10k/".format(data_dir)
		hists=['GH3K27ac_class_data.fs','DH3K27ac_class_data.fs','GH3K27me3_class_data.fs','DH3K27me3_class_data.fs','GH3K4me1_class_data.fs','DH3K4me1_class_data.fs','GH3K9me3_class_data.fs','DH3K9me3_class_data.fs']
		fs_files=[hists_fs_dir+x for x in hists]
		dhm.exists_fs_file=fs_files
		outdir="{}/../Example/delta_bounds_modification/".format(Bin)
		data_dir="{}/data/CHIP_Signal_10k/".format(outdir)
		dhm.outdir=data_dir
		outfig="{}/delta_insulation_bounds_modification.pdf".format(outdir)
		print(outfig)
		fig=plt.figure(figsize=(8,3))
		dhm.draw_all_delta_map(fig)
		fig.savefig(outfig)

if __name__=='__main__':
	unittest.main()
