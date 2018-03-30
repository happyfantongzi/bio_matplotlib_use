#! /usr/bin/env python3
import unittest
import time
import argparse
import numpy as np
import re,sys,os
old_settings = np.seterr(all='ignore')
import pandas as pd
pd.set_option('display.precision', 3)
from ComputMatrix import subMatrix
import matplotlib.gridspec as gridspec
Bin = os.path.abspath(os.path.dirname(__file__))
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')
from scipy import sparse

__author__='fanxuning'
__mail__='xuningfan@genome.cn'


class delta_hic:
	def __init__(self,chr,start,end):
		self.chr=chr
		self.start=start
		self.end=end
		self.sam1=''
		self.sam2=''
		self.mat1_path=''
		self.mat2_path=''
		self.resolution=40000
		self.focus_region=''
		self.up_percent=1
		self.down_percent=99
		self.focusRegion=''
		self.exists_minus_mat=''
		self.outdir=''

	def transform_mat_format(self,path):
		mat = pd.read_table(path,header=0,index_col=0)
		mat=sparse.coo_matrix(mat)
		i,j,data=mat.row,mat.col,mat.data 
		ii=[str(x) for x in i];jj=[str(x) for x in j]
		bin_pairs=pd.Series(ii).str.cat(pd.Series(jj),sep='-')
		transform_mat=pd.concat([pd.Series(i),pd.Series(j),bin_pairs,pd.Series(data)],axis=1)
		transform_mat.columns=['b1','b2','bin_pair','counts']
		transform_mat['counts'].ix[transform_mat['b1']==transform_mat['b2']]=0
		transform_mat=transform_mat.reset_index(drop=True)
		print('transform_mat......')
		print(transform_mat.head())
		return transform_mat
	
	def quantile_normalization(self,df):
		norm_df=df
		sorted_df=pd.DataFrame()
		for col in df:
			sorted_col=pd.Series(sorted(df[col]))
			if sorted_df.shape[0]==0:
				sorted_df=sorted_col
			else:
				sorted_df=pd.concat([sorted_col,sorted_df],axis=1)
		rank = sorted_df.mean(axis = 1).tolist() # 按行取均值
		for col in df:
			t = np.searchsorted(np.sort(df[col]), df[col])
			norm_df[col] = [rank[i] for i in t]
		return norm_df
	
	@property
	def custom_cmap(self):
		my_colors=['#E7F313','#53542A','#040E0F','#007974','#00B8C0']
		custom_cmap = mpl.colors.LinearSegmentedColormap.from_list('cmap',my_colors[::-1],500)
		return custom_cmap

	def minus_matrix(self):
		spearse_mat=pd.DataFrame()
		if not self.exists_minus_mat=='':
			 spearse_mat=np.load(self.exists_minus_mat)
		else:
			control=self.transform_mat_format(self.mat1_path)
			treat=self.transform_mat_format(self.mat2_path)
			df=pd.merge(control[['bin_pair','counts']],treat[['bin_pair','counts']],on='bin_pair',how='outer')
			df=df.reset_index(drop=True)
			df=df.fillna(0)
			df['b1']=df['bin_pair'].str.split('-').str.get(0)
			df['b1']=df['b1'].apply(lambda x: int(x))
			df['b2']=df['bin_pair'].str.split('-').str.get(1)
			df['b2']=df['b2'].apply(lambda x:int(x))
			df.columns=['bin_pair',self.sam1,self.sam2,'b1','b2']
			norm_df=self.quantile_normalization(df[[self.sam1,self.sam2]])
			minus_mat=norm_df[self.sam1]-norm_df[self.sam2]
			condense_mat=pd.concat([df['b1'],df['b2'],minus_mat],axis=1)
			condense_mat.columns=['b1','b2','minus']
			condense_mat['b1']=condense_mat['b1'].apply(lambda x: int(x))
			condense_mat['b2']=condense_mat['b2'].apply(lambda x: int(x))
			out_condense_path="{}/{}_quantile_norm_condense.mat".format(self.outdir,self.chr)
			out_spearse_path="{}/{}_quantile_norm_spearse.npy".format(self.outdir,self.chr)
			condense_mat.to_csv(out_condense_path,header=True,index=False,sep='\t')
			spearse_mat=sparse.coo_matrix((minus_mat,(df['b1'], df['b2'])), dtype=float).toarray()
			np.save(out_spearse_path,spearse_mat)
		spearse_mat=pd.DataFrame(spearse_mat)
		return spearse_mat

	def draw_Interaction_Change_distance_distr(self,ax,percent,spearse_mat):
		strengthen=np.where(spearse_mat>0)
		weaken=np.where(spearse_mat<0)
		strengthen_dist=list(np.abs(strengthen[0]-strengthen[1]))
		weaken_dist=list(np.abs(weaken[0]-weaken[1]))
		ax.tick_params(left='off',right='off',top='off',bottom='off')
		ax.set_title(self.chr)
		labels=['strengthen','weaken']
		bplot=ax.boxplot([strengthen_dist,weaken_dist],
				showfliers=False,labels=labels,patch_artist=True, widths=0.8,vert=True)
		colors=['black', 'white']
		for loc in ['left','bottom']:
			ax.spines[loc].set_linewidth(0.5)
			ax.spines[loc].set_color('k')
		for patch, color in zip(bplot['boxes'], colors):
			patch.set_facecolor(color)
		ax.set_ylabel("distance")

	def draw_minus_matrix(self,ax,raw_mat):
		subM=subMatrix(self.chr,self.start,self.end,self.resolution,self.focusRegion)
		subM.mat=raw_mat
		subM.need_to_be_cut=True
		print('matrix real start , real end:{}\t{}'.format(subM.genome_start,subM.genome_end))
		subM.colormap=self.custom_cmap
		subM.up_percent=self.up_percent
		subM.down_percent=self.down_percent
		#subM.up_down_together=True
		subM.tri_map(ax,subM.all_mat,False)
	
	def draw_heatmap_compared(self,fig):
		#------------  mat1 --------------------
		ax1 = fig.add_axes([0.03,0.07,0.25,0.85],axisbg='w')
		c_subM=subMatrix(self.chr,self.start,self.end,self.resolution,self.focusRegion)
		c_subM.colormap=plt.get_cmap('YlOrRd')
		c_subM.mat_path=self.mat1_path;c_subM.need_to_be_cut=True
		c_subM.sample=self.sam1;c_raw_mat=c_subM.cut_mat(c_subM.all_mat)
		c_subM.plotHeatmap(c_raw_mat,ax1)

		#-------------  mat2 -------------------
		ax2 = fig.add_axes([0.34,0.07,0.25,0.85],axisbg='w')
		t_subM=subMatrix(self.chr,self.start,self.end,self.resolution,self.focusRegion)
		t_subM.mat_path=self.mat2_path;t_subM.need_to_be_cut=True
		t_subM.colormap=plt.get_cmap('YlOrRd')
		t_subM.sample=self.sam2
		t_raw_mat=t_subM.cut_mat(t_subM.all_mat)
		t_subM.plotHeatmap(t_raw_mat,ax2)

		#---------------minus -----------------
		ax3 = fig.add_axes([0.65,0.07,0.25,0.85],axisbg='w')
		m_subM=subMatrix(self.chr,self.start,self.end,self.resolution,self.focusRegion)
		m_subM.mat=self.minus_matrix();m_subM.need_to_be_cut=True
		m_raw_mat=m_subM.cut_mat(m_subM.all_mat)
		m_subM.colormap=self.custom_cmap
		m_subM.sample='{} minus {}'.format(self.sam2,self.sam1)
		m_subM.plotHeatmap(m_raw_mat,ax3)

class TestMatrix(unittest.TestCase):
	def test_delta_hic(self):
		outdir='{}/../Example/delta_Hic/'.format(Bin)
		data_dir="{}/../Example/delta_Hic/data/".format(Bin)
		control_path="{}/GM.mm9_chr19.mm9_chr19.matrix".format(data_dir)
		treat_path='{}/DM.mm9_chr19.mm9_chr19.matrix'.format(data_dir)
		hic_minus=delta_hic('chr19',13000000,33342430)
		hic_minus.sam1='GM'
		hic_minus.sam2='DM'
		hic_minus.mat1_path=control_path
		hic_minus.mat2_path=treat_path
		#hic_minus.exists_minus_mat="/annoroad/data1/bioinfo/PROJECT/Commercial/Cooperation/Hic/ANBJ161163/BJ161163-03/fanxuning/Analysis/ReAnalysis/Analysis_20180312/Diff_boundaries_loc/bin2/Delta_Hic_modification/bin/../Example/delta_Hic//chr19_quantile_norm_spearse.npy"
		hic_minus.outdir=outdir
		'''
		fig=plt.figure(figsize=(16,10))
		ax=fig.add_axes([0.05,0.85,0.8,0.05],axisbg='w')
		raw_mat=hic_minus.minus_matrix()
		hic_minus.draw_minus_matrix(ax,raw_mat)
		fig.savefig('{}/minus_chr19_tri_map.pdf'.format(outdir))
		
		out_heatmap='{}/minus_chr19_heatmap.pdf'.format(outdir)
		print(out_heatmap)
		fig2=plt.figure(figsize=(15,4))
		hic_minus.draw_heatmap_compared(fig2)
		fig2.savefig(out_heatmap)
		'''
		out_boxplot='{}/minus_chr19_boxplot.pdf'.format(outdir)
		print(out_boxplot)
		fig=plt.figure(figsize=(4,6))
		ax=fig.add_axes([0.2,0.1,0.7,0.8],axisbg='white')
		raw_mat=hic_minus.minus_matrix()
		hic_minus.draw_Interaction_Change_distance_distr(ax,0,raw_mat)
		fig.savefig(out_boxplot)

if __name__ == '__main__':
	unittest.main()
