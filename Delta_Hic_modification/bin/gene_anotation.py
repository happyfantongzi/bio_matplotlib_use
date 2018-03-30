#!/usr/bin/env python
import re,sys
import time
import tempfile
import os.path
import sys
import pandas as pd
import numpy as np
Bin = os.path.abspath(os.path.dirname(__file__))
import unittest
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')


__author__='fanxuning'
__mail__='xuningfan@genome.cn'

class Gene_annotation:
	def __init__(self,gtf,chr,start,end):
		self.gtf=gtf
		self.chr=chr
		self.start=start
		self.end=end
		self.show='gene'#transcript
		self.gtf_bed=''
		self.upStream=5000
		self.downStream=5000
		self.name=''
		self.scale=True

	def make_regions(self):
		gtf_bed=pd.read_table(self.gtf,usecols=[0,2,3,4,8],names=['chr','transcript_ID','start','end','desc'])
		gtf_bed=gtf_bed[(gtf_bed['chr']==self.chr) & (gtf_bed['transcript_ID']=='exon')]
		gtf_bed=gtf_bed.reset_index(drop=True)
		gtf_bed=gtf_bed[(gtf_bed['start']>self.start) & (gtf_bed['end']<self.end)]
		gtf_bed=gtf_bed.reset_index(drop=True)
		def gene_match(geneinfo):
			gene_id='';gene_name="";transcript_id='';transcript_name=''
			match=re.search("gene_id \"(.+?)\";.+?gene_name \"(.+?)\";",geneinfo)
			if match:
				gene_id=match.group(1)
				gene_name=match.group(2)
				return (gene_id,gene_name)
	
		def transcript_match(geneinfo):
			transcript_id='';transcript_name=''
			match=re.search("transcript_id \"(.+?)\";.+?transcript_name \"(.+?)\"",geneinfo)
			if match:
				transcript_id=match.group(1)
				transcript_name=match.group(2)
			return (transcript_id,transcript_name)

		gtf_bed['gene_desc']=gtf_bed['desc'].apply(gene_match)
		gtf_bed['transcript_desc']=gtf_bed['desc'].apply(transcript_match)
		self.gtf_bed=gtf_bed

	def draw_genes(self,ax):
		ax.tick_params(left='off',top='off',bottom='off',right='off')
		import matplotlib.patches as mpatches
		for loc in ['top','bottom','left','right']:
			ax.spines[loc].set_visible(False)
		for name,genes in self.gtf_bed.groupby(['gene_desc']): 
			exons=np.array(genes[['start','end']])
			gene_start=np.min(exons.ravel())
			gene_end=np.max(exons.ravel())
			ax.plot([gene_start,gene_end],[0.45,0.45],linewidth=1,color='#143084')
			for exon in exons:
				rec=mpatches.Rectangle((exon[0],0.3),(exon[1]-exon[0]),0.3,color="#143084")
				ax.add_patch(rec)
			#gene_rec=mpatches.Rectangle((gene_start,0.4),(gene_end-gene_start),0.4,ec='#6FCEA3',fc='none')
			#ax.add_patch(gene_rec)
		xmin=self.start-self.upStream;xmax=self.end+self.downStream
		if self.scale:
			ax.annotate('',xy=(self.start,0.98),xytext=(self.end,0.98),fontsize=2,arrowprops=dict(arrowstyle="|-|",color='black',ec='black',linewidth=0.5,connectionstyle="arc3"))
			ax.text((self.start-100),0.98,"{:,}".format(self.start),fontsize=5,horizontalalignment='right',verticalalignment='center',rotation='horizontal',color="black")
			ax.text(0.5,0.8,self.name,transform=ax.transAxes,color='black',fontsize=5)
			ax.text(self.end+100,0.92,"{:,}".format(self.end),fontsize=5,horizontalalignment='left',verticalalignment='center',rotation='horizontal',color="black")
		ax.set_xlim(xmin,xmax)
		xtick=np.linspace(xmin,xmax,5)
		xtlabs=["{:,}".format(x) for x in xtick]
		ax.set_ylim([0,1])
		ax.set_xticks([])
		ax.set_yticks([])

	def draw_transcripts(self,ax):
		import matplotlib.patches as mpatches
		for loc in ['left','right','top','bottom']:
			ax.spines[loc].set_linewidth(1.5)
			ax.spines[loc].set_color('k')
		for trans_id,trans in self.gtf_bed.groupby(['transcript_desc']):
			exons=np.array(trans[['start','end']])
			trans_start=np.min(exons.ravel())
			trans_end=np.max(exons.ravel())
			ax.plot([trans_start,trans_end],[0.6,0.6],linewidth=1,color='#143084')
			for exon in exons:
				rec=mpatches.Rectangle((exon[0],0.4),(exon[1]-exon[0]),0.4,color="#143084")
				ax.add_patch(rec)
			ax.set_ylim([0,1])
			ax.set_xlim(self.start,self.end)
			ax.set_xticks([])
			ax.set_yticks([])

class TestGene_annotation(unittest.TestCase):
	def test_annotation(self):
		data_dir="{}/../Example/gene_annotation/data/".format(Bin)
		outdir="{}/../Example/gene_annotation/".format(Bin)
		print('outfig......')
		outfile="{}/gene_annotation.pdf".format(outdir)
		print(outfile)
		gtfFiles="{}/Mus_musculus.NCBIM37.67.gtf".format(data_dir)
		chr,start,end='chr19',26680000,26720000
		gene_regions=Gene_annotation(gtfFiles,chr,start,end)
		gene_regions.name='myog'
		gene_regions.make_regions()
		fig=plt.figure(figsize=(30,2))
		ax=fig.add_axes([0.1,0.1,0.8,0.4],axisbg='white')
		gene_regions.draw_genes(ax)
		fig.savefig(outfile)

if __name__ == '__main__':
	unittest.main()

