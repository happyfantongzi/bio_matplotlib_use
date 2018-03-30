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
		chr,start,end='chr19',6062821,6068373
		ghm=Gene_modification(gtfFiles,chr,start,end)
		listbams="{}/bamfile.list".format(data_dir)
		ghm.bams_list=listbams
		ghm.name='Tm7sf2'
		hists_fs_dir='/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/bin2/Delta_Hic_modification/bin/../Example/multi_modification/data/multiple_histones_fs_50b/'
		hists=['GH3K27ac_class_data.fs','DH3K27ac_class_data.fs','GH3K27me3_class_data.fs','DH3K27me3_class_data.fs','GH3K4me1_class_data.fs','DH3K4me1_class_data.fs']
		fs_files=[hists_fs_dir+x for x in hists]
		ghm.exists_fs_file=fs_files
		outdir="{}/../Example/gene_modification/".format(Bin)
		data_dir="{}/data/CHIP_Signal_10bp/".format(outdir)
		ghm.outdir=data_dir
		outfig="{}/Afap1l2_modification.pdf".format(outdir)
		print(outfig)
		fig=plt.figure(figsize=(4,4))
		ghm.pre_gene_modication_lay_out()
		ghm.draw_gene_modification(fig)
		fig.savefig(outfig)

if __name__=='__main__':
	unittest.main()
