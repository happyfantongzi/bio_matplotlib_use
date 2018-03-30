#! /usr/bin/env python3
import unittest
import re,sys,os
import numpy as np
import pandas as pd
import matplotlib as mpl
import random
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')
Bin= os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.append(Bin+'/package')
import ZODB
from ZODB import FileStorage, DB
import transaction
from persistent import Persistent
from single_histone_TF import histone_TF
from multiprocessing import Pool
import matplotlib.gridspec as gridspec

__author__='fanxuning'
__mail__='xuningfan@genome.cn'

'''
	this script is used to show multi modification information
	the core is multi processing....
	what's we need?
		list_bams : histone/TF names     treat_bam  control_bam
		genome region: chr,start,end
		outputdir
'''

class MyZODB(Persistent):
	def __init__(self,path):
		self.storage = FileStorage.FileStorage(path)
		self.db = DB(self.storage)
		self.connection = self.db.open()
		self.dbroot = self.connection.root()
	def close(self):
		self.connection.close()
		self.db.close()
		self.storage.close()

class multi_modification:
	def __init__(self):
		self.title=''
		self.chr='chr1'
		self.mstart=0
		self.mend=197195432
		self.bams_list=''
		self.outdir=''
		self.numberOfProcessors=2
		self.data_of_histone_or_TFs=[]
		self.colors=['#953131','#0C1B6C','#64A064','#966161','#914D1D','#7426A4','#75ABFE','#FEC369']
		self.exists_fs_file=None
		self.other_parts_to_show=0
		self.binSize=50
		self.numsOf=0
		self.gs=gridspec.GridSpec(2,1)
		self.modification_layout=[]

	@property
	def listbams(self):
		listbams=pd.read_table(self.bams_list,names=['names','treat','control'])
		return listbams

	def pre_modification_lay_out(self):
		print('pre lay out.......')
		'''
			use gridspect to set lay out
		'''
		other_parts=self.other_parts_to_show
		subplots=[];numsOf=self.listbams.shape[0]+other_parts
		gs = gridspec.GridSpec(numsOf, 1)
		gs.update(hspace=0.7)
		for i in range(other_parts,numsOf):
			self.modification_layout.append(i)
		self.numsOf=numsOf
		self.gs=gs

	def multi_processing(self):
		listbams=self.listbams
		pool = Pool(self.numberOfProcessors)
		for i,r in listbams.iterrows():
			outputfile="{}/{}_log2ratio.bed".format(self.outdir,r['names'])
			protein=histone_TF(r['names'])
			protein.chr=self.chr
			protein.start=self.mstart
			protein.end=self.mend
			protein.treat_bam=r['treat']
			protein.control_bam=r['control']
			protein.binSize=self.binSize
			protein.outFileName=outputfile
			pool.map_async(protein.deeptools_bamCompare(),[]).get()
			protein_dbfile="{}/{}_class_data.fs".format(self.outdir,r['names'])
			protein.outdbfile=protein_dbfile
			self.data_of_histone_or_TFs.append(protein_dbfile)
			db = MyZODB(protein_dbfile)
			dbroot=db.dbroot
			dbroot['load']=protein
			transaction.commit()
			db.close()
		pool.close() 
		pool.join()

	def multi_draw_ChIP_signal(self,fig):
		print('draw multi chip sigal.....')
		#subplots=self.pre_lay_out()
		if self.exists_fs_file:
			self.data_of_histone_or_TFs=self.exists_fs_file
			print('len....')
			print(len(self.exists_fs_file))
		else:
			self.multi_processing()
		print('self.modification_layout.....')
		print(self.modification_layout)
		for i,prodb in enumerate(self.data_of_histone_or_TFs):
			color=random.choice(self.colors)
			db = MyZODB(prodb)
			dbroot=db.dbroot
			pro=dbroot['load']
			pro.start=self.mstart
			pro.end=self.mend
			gridsSub= gridspec.GridSpecFromSubplotSpec(1,1,
					subplot_spec=self.gs[self.modification_layout[i],:],wspace=0.5,hspace=0.7)
			ax=fig.add_subplot(gridsSub[:,:],axisbg='white')
			#ax = plt.subplot(gs[self.modification_layout[i], :],axisbg='white')
			barplot=pro.draw_modification_coverage(ax,color)
			transaction.commit()
			db.close()

	def multi_peak_draw(self):
		print('draw multi peaks.......')
		#TODO

class multi_modification_test(unittest.TestCase):
	def test_multi_modification(self):
		listbams="{}/../Example/multi_modification/data/bamfile.list".format(Bin)
		outdir="{}/../Example/multi_modification/".format(Bin)
		#testcase1
		test_dir="{}/../Example/multi_modification/data/multiple_histones_fs_50b/".format(Bin)
		print('out_fs_dir.....')
		print(test_dir)
		hists=['H3K27ac_class_data.fs','H3K27me3_class_data.fs','H3K4me1_class_data.fs','H3K9me3_class_data.fs']
		fs_files=[test_dir+x for x in hists]
		mmd=multi_modification()
		mmd.chr='chr19'
		mmd.start=33000000
		mmd.end=53342430
		mmd.bams_list=listbams
		mmd.outdir=test_dir
		#mmd.exists_fs_file=fs_files
		#mmd.other_parts_to_show=5
		#testcase2
		fig=plt.figure(figsize=(4,1))
		mmd.pre_modification_lay_out()
		mmd.multi_draw_ChIP_signal(fig)
		out_fig="{}/multi_chip_signal.pdf".format(outdir)
		print(out_fig)
		fig.savefig(out_fig)

if __name__=='__main__':
	unittest.main()
