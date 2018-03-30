#!/usr/bin/env python3
import re,sys,os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')
Bin= os.path.dirname(os.path.abspath(sys.argv[0]))
import unittest
from gene_anotation import Gene_annotation
from multi_modification import multi_modification
import matplotlib.gridspec as gridspec

__author__='fanxuning'
__mail__='xuningfan@genome.cn'

class Gene_modification(Gene_annotation,multi_modification):
	def __init__(self,gtf,chr,start,end):
		Gene_annotation.__init__(self,gtf,chr,start,end)
		multi_modification.__init__(self)
		self.gene_grids=1
		self.mstart=self.start+self.upStream
		self.mend=self.end+self.downStream

	def pre_gene_modication_lay_out(self):
		'''
			1.标记基因区域大小
			2.画外显子
			gene grid  : 1
			multi modification : 8
		'''
		self.pre_modification_lay_out()
		all_grids=self.gene_grids+self.listbams.shape[0]+1
		Chip_layouts=list(range(self.gene_grids,all_grids))
		self.gs=gridspec.GridSpec(all_grids,1)
		self.modification_layout=Chip_layouts

	def draw_gene_modification(self,fig):
		#--------------------------draw  genes -----------------------------------------#
		G=self.gs
		gridsSub1=gridspec.GridSpecFromSubplotSpec(1,1,
			subplot_spec=G[0:self.gene_grids,:],wspace=0.5,hspace=0.7)
		ax1=fig.add_subplot(gridsSub1[:,:],axisbg='white')
		self.make_regions()
		self.draw_genes(ax1)
		#--------------------------draw histone  modification---------------------------#
		self.multi_draw_ChIP_signal(fig)

class Test_Gene_modification(unittest.TestCase):
	def test_gene_modification(self):
		data_dir="{}/../Example/gene_modification/data/".format(Bin)
		gtfFiles="{}/Mus_musculus.NCBIM37.67.gtf".format(data_dir)
		chr,start,end='chr19',56948739,56986568
		ghm=Gene_modification(gtfFiles,chr,start,end)
		listbams="{}/bamfile.list".format(data_dir)
		ghm.bams_list=listbams
		ghm.name='Vwa2'
		hists_fs_dir="{}/multiple_histones_fs_5b/".format(data_dir)
		hists=['GH3K27ac_class_data.fs','DH3K27ac_class_data.fs','GH3K27me3_class_data.fs','DH3K27me3_class_data.fs','GH3K4me1_class_data.fs','DH3K4me1_class_data.fs']
		fs_files=[hists_fs_dir+x for x in hists]
		ghm.exists_fs_file=fs_files
		outdir="{}/../Example/Delta_Hic_modification/".format(Bin)
		chip_data_dir="{}/data/CHIP_Signal_5bp/".format(outdir)
		ghm.outdir=chip_data_dir
		outfig="{}/Vwa2_modification.pdf".format(outdir)
		fig=plt.figure(figsize=(4,4))
		ghm.pre_gene_modication_lay_out()
		ghm.draw_gene_modification(fig)
		fig.savefig(outfig)

if __name__=='__main__':
	unittest.main()
