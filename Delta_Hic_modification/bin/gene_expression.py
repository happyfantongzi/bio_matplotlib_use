#! /usr/bin/env python3
import os,re,sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import unittest
from matplotlib import pyplot as plt
plt.style.use('ggplot')
Bin = os.path.abspath(os.path.dirname(__file__))

class Gene_Express:
	def __init__(self,chr,start,end):
		self.chr=chr
		self.start=start
		self.end=end
		self.annote_file=""

	@property
	def extract_genes(self):
		anno_df = pd.read_table(self.annote_file,header=0,index_col=None,encoding='utf-8')
		position_df=anno_df['Position'].str.split(":|-",n=3,expand=True)
		position_df.columns=['chr','start','end','strand']
		anno_df=pd.concat([anno_df,position_df],axis=1)
		anno_df['start']=anno_df['start'].apply(lambda x:int(x))
		anno_df['end']=anno_df['end'].apply(lambda x: int(x))
		sub_anno_df=anno_df[(anno_df['chr']==self.chr)&(anno_df['start']>self.start-1)&(anno_df['end']<self.end+1)]
		sub_anno_df=sub_anno_df.reset_index(drop=True)
		sub_anno_df=sub_anno_df[['chr','start','end','Gene','FoldChange','Significant']]
		sub_anno_df=sub_anno_df.sort_values(by='start')
		#sub_anno_df=sub_anno_df[sub_anno_df['Significant']=='yes']
		sub_anno_df=sub_anno_df.reset_index(drop=True)
		sub_anno_df=sub_anno_df.replace(np.inf, 10)
		sub_anno_df=sub_anno_df.fillna(0)
		sub_anno_df=sub_anno_df.replace(-np.inf,-10)
		return sub_anno_df
	
	def draw_gene_expression(self,ax):
		ax.tick_params(left='off',top='off',bottom='on',right='off')
		for loc in ['top','bottom','left','right']:
			ax.spines[loc].set_linewidth(0.5)
			ax.spines[loc].set_color('k')
		#ax.spines['bottom'].set_position(('outward',2))
		for i,r in self.extract_genes.iterrows():
			width=r['end']-r['start']
			color='red' if r['FoldChange']>0 else 'blue'
			ax.bar(left=r['start'], height=r['FoldChange'], width=width, 
					bottom=0,color=color,align="edge",edgecolor=color)
		ax.text(-0.01,0.5,'GeneExpress',color='k',horizontalalignment='right',verticalalignment='center',rotation='horizontal',transform=ax.transAxes,fontsize=4)
		ax.set_ylim(-10,10)
		xticks=list(np.linspace(self.start,self.end,6))
		print('start')
		print(self.start)
		print('end')
		print(self.end)
		xticklabs=list(np.linspace(self.start,self.end,6))
		tkls=['{:,}'.format(int(i)) for i in xticklabs]
		print('xticks....')
		print(xticks)
		print(tkls)
		ax.set_xticks(xticks)
		ax.set_xticklabels(tkls,fontsize=4)
		ax.set_xlim(self.start,self.end)
		ax.set_yticks([])
		

class test_Gene_express(unittest.TestCase):
	def test_gene_express(self):
		gene_expr=Gene_Express('chr19',10000000,48342430)
		outdir="{}/../Example/Gene_express/".format(Bin)
		annote_file="{}/data/DM_GM.anno.xls".format(outdir)
		gene_expr.annote_file=annote_file
		outfig="{}/gene_expression.pdf".format(outdir)
		print(outfig)
		fig=plt.figure(figsize=(10,1))
		ax=fig.add_axes([0.1,0.1,0.8,0.8],axisbg='w')
		gene_expr.draw_gene_expression(ax)
		fig.savefig(outfig)

if __name__=='__main__':
	unittest.main()
