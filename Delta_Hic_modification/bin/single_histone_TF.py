#! /usr/bin/env python3
import argparse
import pyBigWig
import unittest
import re,sys,os
import numpy as np
import pandas as pd
import  matplotlib 
matplotlib.use('Agg')
import random
Bin= os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.append(Bin+'/package')
from persistent import Persistent
from matplotlib import pyplot as plt
plt.style.use('ggplot')
plt.rcParams['font.size']='12'
old_settings = np.seterr(all='ignore')
'''
	this wrap up deeptools to show  histone/TF  modification 
	what's need ?
			sample name 
			histone/TF name
			bamfiles 
			peakfiles
			binsize 
'''

class histone_TF(Persistent):
	def __init__(self,name):
		self.name=name
		self.treat_bam=''
		self.control_bam=''
		self.peak=''
		self.binSize=40000
		self.chr='chr1'
		self.start=0
		self.end=197195432
		self.focus_region=''
		self.justCoverage=False
		self.needCompare=True
		self.outFileName=''
		self.outdbfile=''
		self.numberOfProcessors=2
		self.show_focus_regions=False

	def deeptools_bamCoverage(self):
		func_args = {'scaleFactor': 1.0}
		from deeptools import bamCoverage
		from deeptools import writeBedGraph
		wr=writeBedGraph.WriteBedGraph([self.treat_bam,self.control_bam],
				binLength=self.binSize,
				stepSize=self.binSize,
				region=None,
				blackListFileName=None,
				numberOfProcessors=self.numberOfProcessors,
				extendReads=False,
				minMappingQuality=None,
				ignoreDuplicates=True,
				center_read=False,
				zerosToNans=False,
				samFlag_include=None,
				samFlag_exclude=None,
				minFragmentLength=0,
				maxFragmentLength=0,
				verbose=True,
				)
		wr.run(writeBedGraph.scaleCoverage, func_args, self.outFileName,
				blackListFileName=None,
				)

	def deeptools_bamCompare(self):
		print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
		print(self.name)
		from deeptools import writeBedGraph
		from deeptools.SES_scaleFactor import estimateScaleFactor
		from deeptools import parserCommon
		from deeptools import bamHandler
		from deeptools import getRatio
		from deeptools.getScaleFactor import get_num_kept_reads
		debug = 0
		FUNC = getRatio.getRatio
		#log2,ratio,subtract,add,mean,reciprocal_ratio,first,second
		func_args = {'valueType': "ratio",
				'scaleFactors': [1,1],
				'pseudocount': 1
				}
		from deeptools import writeBedGraph
		wr = writeBedGraph.WriteBedGraph([self.treat_bam,self.control_bam],self.binSize,0,
				stepSize=self.binSize,
				region=None,
				numberOfProcessors=self.numberOfProcessors,
				extendReads=False,
				blackListFileName=None,
				minMappingQuality=None,
				ignoreDuplicates=False,
				center_read=False,
				zerosToNans=False,
				samFlag_include=None,
				samFlag_exclude=None,
				minFragmentLength=0,
				maxFragmentLength=0,
				verbose=True,
				)
		wr.run(FUNC,func_args,self.outFileName)

	def make_regions(self):
		print("++++++++++++++++++++++++++++++++++ outfile +++++++++++++++++++++++++++++++++++++")
		print(self.outFileName)
		reads_counts=pd.read_table(self.outFileName,names=['chr','start','end','signal'])
		reads_counts=reads_counts[reads_counts['chr']==self.chr]
		reads_counts=reads_counts.reset_index(drop=True)
		reads_counts=reads_counts[(reads_counts['start']>(self.start-1))&(reads_counts['end']<(self.end+1))]
		reads_counts=reads_counts.reset_index(drop=True)
		return reads_counts

	def draw_modification_coverage(self,ax,color):
		reads_counts=self.make_regions()
		max=np.max(np.array(reads_counts['signal']))
		min=np.min(np.array(reads_counts['signal']))
		ax.tick_params(left='off',top='off',right='off',bottom='off')
		#----------set spines-------------------------------------------#
		for loc in ['right','top','bottom']:
			ax.spines[loc].set_visible(False)
		ax.spines['left'].set_linewidth(0.5)
		#ax.spines['left'].set_position(('outward',30))
		ax.spines['left'].set_color('k')
		#---------bar  plot --------------------------------------------#
		barplot=ax.bar(left=reads_counts['start'], height=reads_counts['signal'], width=(self.binSize)*0.8, bottom=0,color=color,align="edge",edgecolor=color)
		#---------set labels -------------------------------------------#
		ax.set_xticks([])
		ax.set_ylim(1,max)
		ax.set_yticks([])
		'''
		ax.set_yticks([1,max])
		ax.set_yticklabels([1,max],color='red',fontsize=5)
		'''
		print('{} start,end:{}\t{}'.format(self.name,self.start,self.end))
		ax.set_xlim(self.start,self.end)
		ax.text(-0.01,0.5,self.name,horizontalalignment='right',verticalalignment='center',transform=ax.transAxes,fontsize=4)
		#--------display precision-------------#
		import matplotlib.ticker as ticker
		ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
		return barplot

	def draw_peaks(self):
		print('draw peaks.....')
		listpeak=pd.read_table(self.peak_list,names=['histone','path'])
		for i,r in listpeak.iterrows(): 
			df=pd.read_table(r['path'],names=['chrom','start','end','name','fold'])
			df=df[(df['start']>(self.real_start-1)) & (df['end']<(self.real_end+1))]
			df=df.reset_index()
			describe = df['fold'].describe()
			df['fold']=(df['fold']/describe['max'])*0.7
			ax.text(self.real_start,i,r['histone'],color='k',horizontalalignment='right',verticalalignment='center',fontsize=2)
			plot=ax.bar(left=df['start'], height=df['reads'], width=(self.resolution)*0.8, bottom=i,color=color,align="edge",edgecolor=color)
			ax.text(-0.05,0.5,'CHIP seq',color='k',horizontalalignment='right',verticalalignment='center',rotation='vertical',transform=ax.transAxes)


class Testmodification(unittest.TestCase):
	def test_modification(self):
		colors=['#953131','#0C1B6C','#64A064','#966161','#914D1D','#7426A4','#75ABFE','#FEC369']
		color=random.choice(colors)
		name='H3K27ac'
		treat='/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/EPI/Chip/ANBJ161163/BJ161163-10/jintiantian/Analysis/Analysis/Analysis/Alignment/D_K27ac/D_K27ac_1/D_K27ac_1.sort.bam'
		control='/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/EPI/Chip/ANBJ161163/BJ161163-10/jintiantian/Analysis/Analysis/Analysis/Alignment/DM_INPUT_4/DM_INPUT_4/DM_INPUT_4.sort.bam'
		'''
		treat='{}/../Example/single_histone_TF/data/D_K27ac_1.sort.bam'.format(Bin)
		control='{}/../Example/single_histone_TF/data/DM_INPUT_4.sort.bam'.format(Bin)
		'''
		peak='{}/../Example/single_histone_TF/data/D_K27ac_peaks.bed'.format(Bin)
		H3K27ac=histone_TF(name)
		H3K27ac.control_bam=control
		H3K27ac.treat_bam=treat
		#H3K27ac.binSize=10000
		H3K27ac.outFileName="{}/../Example/single_histone_TF/H3K27ac_log2ratio.bed".format(Bin)
		outfig='{}/../Example/single_histone_TF/H3K27ac_Chip_signal.pdf'.format(Bin)
		#H3K27ac.deeptools_bamCoverage()
		H3K27ac.deeptools_bamCompare()
		fig=plt.figure(figsize=(16,1))
		ax=fig.add_axes([0.2,0.1,0.8,0.8],axisbg='white')
		H3K27ac.draw_modification_coverage(ax,color)
		fig.savefig(outfig)

if __name__ == '__main__':
	unittest.main()

