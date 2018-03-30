'''
	a simple class  store sample info 
	need:
		genomic region : chr start end
		extend : beforeRegionLength  afterRegionStartLength
		normlized matrix : ice matrix better 
		TAD info :  insulation  boundaries   
		gene info : bed file better
		peak info : bw file list or peak list,peak  better
		outdir : output dir 

	notice : need to check out start and end
'''
import re,sys,os,math 
Bin = os.path.abspath(os.path.dirname(__file__))
import numpy as np
import pandas as pd
import unittest
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')
import matplotlib.gridspec as gridspec
from matplotlib import colors
from single_insulation_TAD import TAD_insulation
from multi_modification import multi_modification

class Sample(TAD_insulation,multi_modification):

	def __init__(self,name,chr,start,end):
		TAD_insulation.__init__(self,name,chr,start,end)
		multi_modification.__init__(self)

	def pre_all_lay_out(self):
		print("numof:\t{}\n".format(self.numsOf))
		print('print all_lay_out.....')
		self.prepare_TAD_lay_out()
		self.pre_modification_lay_out()
		all_grids=self.grids+self.listbams.shape[0]
		Chip_layouts=list(range(self.grids,all_grids))
		self.modification_layout=Chip_layouts
		self.Gs=gridspec.GridSpec(all_grids,1)
		self.gs=gridspec.GridSpec(all_grids,1)

	def draw_sample_info(self,fig):  
		#self.multi_processing()
		self.pre_all_lay_out()
		self.draw_TAD(fig)
		self.mstart=self.real_start
		self.mend=self.real_end
		print('real start ...real end....')
		print(self.mstart,self.mend)
		self.pre_all_lay_out()
		self.multi_draw_ChIP_signal(fig)

class TestSample(unittest.TestCase):
	def test_sample(self):
		sam=Sample('DM','chr19',30000000,45342430)
		listbams="{}/../Example/Sample_info/data/DM_bamfile.list".format(Bin)
		print('listbams....')
		print(listbams)
		sam.mat_path='{}/../Example/Sample_info/data/DM1.mm9_chr19.mm9_chr19.matrix'.format(Bin)
		print('mat_path')
		print(sam.mat_path)
		sam.boundaries='{}/../Example/Sample_info/data/DM1_chr19_immeran.boundaries'.format(Bin)
		print('sam.boundaries......')
		print(sam.boundaries)
		sam.insulation='{}/../Example/Sample_info/data/DM1_chr19_immeran.insulation'.format(Bin)
		print('sam.insulation....')
		print(sam.insulation)
		sam.outdir="{}/../Example/Sample_info/".format(Bin)
		print('outdir....')
		print(sam.outdir)
		out_fig='{}/sample_info.pdf'.format(sam.outdir)

		test_dir="{}/../Example/multi_modification/data/multiple_histones_fs/".format(Bin)
		hists=['H3K27ac_class_data.fs','H3K27me3_class_data.fs','H3K4me1_class_data.fs','H3K9me3_class_data.fs']
		fs_files=[test_dir+x for x in hists]
		sam.exists_fs_file=fs_files
		sam.resolution=40000
		sam.other_parts_to_show=5
		sam.extend=''
		sam.region=''
		sam.bams_list=listbams
		fig=plt.figure(figsize=(5,2))
		#G=sam.pre_all_lay_out()
		sam.draw_sample_info(fig)
		fig.savefig(out_fig)
		print(out_fig)

if __name__ == '__main__':
	unittest.main()
