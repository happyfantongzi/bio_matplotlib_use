#! /usr/bin/env python3
import re,sys,os
import numpy as np
import pandas as pd
import matplotlib 
matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')
from single_insulation_TAD import TAD_insulation
Bin=os.path.dirname(os.path.abspath(sys.argv[0]))
sys.path.append(Bin+'/package')
from classify_tad_bounds import Bounds_classifier
import matplotlib.gridspec as gridspec
import unittest
import matplotlib.ticker as ticker
'''
	this is two samples boundaries minus 
	what's neeed?
		boundary files
		insulation files
'''
class delta_boundaries(): 
	def __init__(self):
		self.sam1=''
		self.sam2=''
		self.bound1_path=''
		self.bound2_path=''
		self.insulat1=''
		self.insulat2=''
		self.chr=''
		self.start=0
		self.end=0
		self.genome_start=0
		self.genome_end=0
		self.interval=80000
		self.bounds_classify=pd.DataFrame()
		self.color1='#372C78'
		self.color2='#CECFD3'
		self.flag=False
		self.bi_Gs=gridspec.GridSpec(10,1)
		self.bi_layout=[]
		self.show_insulations=True
		self.show_delta_insulation=True
		self.show_boundaries1=True
		self.show_boundaries2=True
		self.title=''

	def sample_bounds(self,bound_path):
		sam_bd=TAD_insulation(self.sam1,self.chr,self.start,self.end)
		sam_bd.boundaries=bound_path
		return sam_bd

	def classfied_boundaries(self):
		bc=Bounds_classifier(self.sam1,self.sam2)
		bo1=self.sample_bounds(self.bound1_path)
		bo2=self.sample_bounds(self.bound2_path)
		bc.bound1=bo1.all_bounds
		bc.bound2=bo2.all_bounds
		#this common bounds is a dict and the interval is 80000
		bounds_classify=bc.classify_boundaries()
		print('classfied boundaries real_start,real_end:{}\t{}'.format(bo1.real_start,bo1.real_end))
		self.genome_start=bo1.real_start
		self.genome_end=bo1.real_end
		bounds_classify=bounds_classify[(bounds_classify['start']>bo1.real_start-1)
				& (bounds_classify['end']<bo1.real_end+1)]
		start_end_df=pd.DataFrame([[bo1.real_start,(bo1.real_start+40000),0,'start'],
			[bo1.real_end,bo1.real_end+40000,0,'end']],
			columns=['start','end','boundaryStrength','type'])
		print(start_end_df)
		bounds_classify=bounds_classify.append(start_end_df)
		bounds_classify=bounds_classify.sort_values(by='start')
		print(bounds_classify.head())
		bounds_classify=bounds_classify.reset_index(drop=True)
		self.bounds_classify=bounds_classify

	def draw_classify_domain(self,ax):
		ax.tick_params(left='off',top='on',bottom='off',right='off')
		print('draw diff boundaries....')
		'''
		three kind of boundaries:
			overlaping
			dispearing 虚线标识
			reforming  实线标识
		'''
		import matplotlib.patches as mpatches
		widths=list(np.diff(list(self.bounds_classify['start'])))
		widths.append(0)
		for i,r in self.bounds_classify.iterrows():
			y=0.1 if (i%2==0) else 0.5;x=r['start']
			color=self.color1 if (i%2==0) else self.color2
			rec=mpatches.Rectangle((x,y),widths[i],0.4, color=color)
			if r['type']=='dispearing':
				ax.vlines(r['start'], [0],1,lw=4,color="grey",alpha=0.5)
			if r['type']=='reforming':
				ax.vlines(r['start'],[0],1,lw=4,color="#00FFFF",linestyle='--',alpha=0.5)
			ax.add_patch(rec)
		ax.set_ylim([0,1])
		ax.set_yticks([])
		if self.flag:
			ax.set_xticks(np.linspace(0,(self.end-self.start),5))
			ax.set_xticklabels(np.linspace(self.start,self.end,5))
		else:
			ax.set_xticks([])
		ax.set_xlim(self.start,self.end)
		###notice  ax.set_xlim() must be consistent.......

	def draw_insulations(self,ax):
		print('#---------------- draw insulations ------------------------------ #')
		bo1=self.sample_bounds(self.bound1_path)
		bo1.insulation=self.insulat1;bo1.icolor='red'
		bo2=self.sample_bounds(self.bound2_path)
		bo2.insulation=self.insulat2;bo2.icolor='green'
		bo1.draw_insulation(ax)
		bo2.draw_insulation(ax)
		ax.set_xticks([])
		ax.set_yticks([])
		ax.text(-0.01,0.5,'insulation',horizontalalignment='right',verticalalignment='center',transform=ax.transAxes,fontsize=4)
	
	def minus_insulations(self):
		print('minus insulation........')
		bo1=self.sample_bounds(self.bound1_path)
		bo1.insulation=self.insulat1
		bo2=self.sample_bounds(self.bound2_path)
		bo2.insulation=self.insulat2
		c_ins_df=bo1.cut_insulation();t_ins_df=bo2.cut_insulation()
		minus_df=c_ins_df.merge(t_ins_df,suffixes=['_{}'.format(self.sam1), '_{}'.format(self.sam2)],on='mid')
		minus_df=minus_df.fillna(0)
		minus_df['minus']=minus_df['insulationScore_{}'.format(self.sam2)]-minus_df['insulationScore_{}'.format(self.sam1)]
		print(minus_df.head())
		return minus_df

	def draw_delta_insulation(self,ax):
		minus_df=self.minus_insulations()
		print('draw compared insulation....')
		for loc in ['top','left','right','bottom']:
			ax.spines[loc].set_visible(False)
		ax.tick_params(left='off',right='off',top='off',bottom='off')
		ax.fill_between(minus_df['mid'],0, minus_df['minus'],where=minus_df['minus']>=0,facecolor='black')
		ax.fill_between(minus_df['mid'],0,minus_df['minus'],where=minus_df['minus']<0,facecolor='gray')
		ax.set_xlim(self.start,self.end)
		ymax=max(minus_df['minus'])
		ymin=min(minus_df['minus'])
		ax.set_xticks([])
		ax.set_yticks([ymin,ymax])
		ax.set_yticklabels([ymin,ymax],fontsize=3)
		ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
		ax.set_ylim(ymin*1.5,ymax*1.5)
		ax.text(-0.01,0.5,'delta insulation',horizontalalignment='right',verticalalignment='center',transform=ax.transAxes,fontsize=4)

	def pre_delta_bound_layout(self):
		'''
		layout set ......
		draw insulation : 2 grids
		delta insulation : 2 grids
		boundaries1:  1 grid
		boundaries2:  1 grid
		'''
		bi_layout=[1,]
		grids=1
		if self.show_insulations:
			grids+=2;bi_layout.append(grids)
		if self.show_delta_insulation:
			grids+=2;bi_layout.append(grids)
		if self.show_boundaries1: 
			grids+=1;bi_layout.append(grids)
		if self.show_boundaries2:
			grids+=1;bi_layout.append(grids)
		self.bi_layout=bi_layout[::-1]
		Gs=gridspec.GridSpec(grids,1)
		self.bi_Gs=Gs
		return grids

	def draw_delta_map(self,fig):
		G=self.bi_Gs
		l=self.bi_layout
		plt.suptitle(self.title,fontsize=12)
		#----------------draw insulations ----------------------#
		if self.show_insulations:
			gridsSub = gridspec.GridSpecFromSubplotSpec(1,1, subplot_spec=G[l.pop():l[-1],:], wspace=0.5, hspace=0)
			ax1=fig.add_subplot(gridsSub[:,:],axisbg='white')
			ax1.set_yticks([])
			self.draw_insulations(ax1)

		#-------------draw  delta insulation -------------------#
		if self.show_delta_insulation:
			gridsSub=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=G[l.pop():l[-1]],wspace=0.5,hspace=0)
			ax2=fig.add_subplot(gridsSub[:,:],axisbg='white')
			self.draw_delta_insulation(ax2)

		#------------draw boundary 1 ----------------------------#
		if self.show_boundaries1:
			gridSub=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=G[l.pop():l[-1]],wspace=0.5,hspace=0)
			ax3=fig.add_subplot(gridSub[:,:],axisbg='white')
			bo1=self.sample_bounds(self.bound1_path)
			bo1.draw_bounds(ax3)
			ax3.spines['top'].set_visible(False)
			ax3.text(-0.01,0.5,'{} bounds'.format(self.sam1),horizontalalignment='right',verticalalignment='center',transform=ax3.transAxes,fontsize=4)

		#------------draw boundary 2 ----------------------------#
		if self.show_boundaries2:
			gridSub=gridspec.GridSpecFromSubplotSpec(1,1,subplot_spec=G[l.pop():l[-1]],wspace=0.5,hspace=0)
			ax4=fig.add_subplot(gridSub[:,:],axisbg='white')
			bo2=self.sample_bounds(self.bound2_path)
			bo2.draw_bounds(ax4)
			ax4.spines['top'].set_visible(False)
			ax4.text(-0.01,0.5,"{} bounds".format(self.sam2),horizontalalignment='right',verticalalignment='center',transform=ax4.transAxes,fontsize=4)


class test_delta_bounds(unittest.TestCase):
	def test_delta_bounds(self): 
		data_path='{}/../Example/delta_tad_bounds/data/'.format(Bin)
		outdir="{}/../Example/delta_tad_bounds/".format(Bin)
		bound1_path="{}/GM_chr19.boundaries".format(data_path)
		bound2_path="{}/DM_chr19.boundaries".format(data_path)
		insulat1_path="{}/GM_chr19.insulation".format(data_path)
		insulat2_path="{}/DM_chr19.insulation".format(data_path)
		debs=delta_boundaries()
		debs.bound1_path=bound1_path
		debs.bound2_path=bound2_path
		debs.insulat1=insulat1_path
		debs.insulat2=insulat2_path
		debs.chr='chr19'
		debs.sam1='GM'
		debs.sam2='DM'
		#debs.start=33000000
		debs.start=3300000
		debs.end=59342430
		debs.classfied_boundaries()
		fig=plt.figure(figsize=(10,1))
		#tb.draw_bounds(ax)
		#debs.draw_classify_domain(ax)
		debs.pre_delta_bound_layout()
		debs.draw_delta_map(fig)
		print('outfig....')
		print("{}/delta_tad.pdf".format(outdir))
		fig.savefig("{}/delta_tad.pdf".format(outdir))

if __name__=='__main__':
	unittest.main()
