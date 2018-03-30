'''
	a simple class  store TAD  
	need:
		genomic region : chr start end
		extend : beforeRegionLength  afterRegionStartLength
		normlized matrix : ice matrix better 
		TAD info :  insulation  boundaries   
		gene info : bed file better
		peak info : bw file list or peak list,peak  better
		outdir : output dir 
		update : TAD grids from 3 to 2
'''
from ComputMatrix import subMatrix
import re,sys,os,math 
import numpy as np
import pandas as pd
import unittest
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
plt.style.use('ggplot')
import matplotlib.gridspec as gridspec
from matplotlib import colors
Bin = os.path.abspath(os.path.dirname(__file__))

class TAD_insulation:
	def __init__(self,name,chr,start,end):
		self.name=name
		self.chr=chr
		self.start=start
		self.end=end
		self.extend=''
		self.region=''
		self.mat_path=''
		self.mat=pd.DataFrame()
		self.insulation=''
		self.boundaries=''
		self.resolution=40000
		self.outdir=''
		self.interval=80000
		self.bounds_need_to_be_merged=True
		self.show_insulation=True
		self.show_boundaries=True
		self.colormap=plt.get_cmap('Reds')
		self.show_map=True
		colorslist=['#004A1E','#005D25','#097532','#2B934B','#67BD6F','#85CC84']
		self.bounds_cmap=colors.LinearSegmentedColormap.from_list('mycmap',colorslist[::-1],100)
		self.title="{} TAD dekker".format(name)
		self.grids=0
		self.Gs=gridspec.GridSpec(2,1)
		self.tad_grid_layout=[]
		self.icolor='k'

	@property
	def sample_mat(self):
		mat=subMatrix(self.chr,self.start,self.end,self.resolution,self.region)
		mat.mat=self.mat
		mat.need_to_be_cut=True
		mat.mat_path=self.mat_path
		mat.colormap=self.colormap
		return mat
	
	@property
	def real_start(self):
		start=int(self.start/self.resolution)*self.resolution
		return start

	@property
	def real_end(self):
		end=math.ceil(self.end/self.resolution)*self.resolution
		return end
	
	@property
	def all_bounds(self):
		bound_df=pd.read_table(self.boundaries,header=0,index_col=None,encoding="utf-8")
		bound_df=bound_df.fillna(0)
		bound_df=bound_df.loc[:,['start','end','boundaryStrength']]
		print('boundaries......')
		print(bound_df.head())
		if self.bounds_need_to_be_merged:
			self.merge_nearby_borders(bound_df)
		return bound_df
	
	def merge_nearby_borders(self,bound_df):
		nearby_dist=np.diff(list(bound_df['start']))
		nearby=np.where(nearby_dist<=self.interval)
		new_starts=[list(bound_df.loc[x:x+1,:].mean()) for x in nearby[0]]
		drop_rows=[x+1 for x in nearby[0]];drop_rows.extend(list(nearby[0]))
		bound_df=bound_df.drop(drop_rows)
		bound_df=bound_df.reset_index(drop=True)
		bound_df=bound_df.append(pd.DataFrame(new_starts,columns=bound_df.columns))
		bound_df=bound_df.sort_values(by='start')
		bound_df=bound_df.reset_index(drop=True)

	def cut_bound(self):
		print("TAD insulation start:end\n{}\t{}".format(self.real_start,self.real_end))
		bound_df=self.all_bounds
		bound_df['mid']=(bound_df['start']+bound_df['end'])/2
		bound_df=bound_df[(bound_df['start']>self.real_start-1)&(bound_df['start']<self.real_end+1)]
		bound_df=bound_df.reset_index(drop=True)
		return bound_df

	def draw_bounds(self,ax):
		bound_df=self.cut_bound()
		ax.tick_params(top='off',bottom='off',left='off',right='off')
		#----------------- modified as  scatter------------------------#
		y=np.zeros(bound_df.shape[0])
		vmin=np.min(np.array(bound_df['boundaryStrength']))
		vmax=np.max(np.array(bound_df['boundaryStrength']))
		norm=colors.Normalize(vmin=vmin, vmax=vmax)
		ax.scatter(bound_df['mid'],y,c=bound_df['boundaryStrength'],cmap=self.bounds_cmap,norm=norm,edgecolor='None',s=500,lw=1.5,marker="|")
		ax.set_xticks([])
		ax.set_yticks([])
		ax.set_xlim(self.real_start,self.real_end)
		for loc in ['left','top','right','bottom']:
			ax.spines[loc].set_visible(False)
			#ax.spines[loc].set_linewidth(1.5)
			#ax.spines[loc].set_color('k')

	def draw_insulation(self,ax):
		import matplotlib.ticker as ticker
		ins_mid=self.cut_insulation()
		ax.tick_params(top='off',bottom='off',left='on',right='off')
		line=ax.plot(ins_mid['mid'], ins_mid['insulationScore'], color=self.icolor, linewidth=0.5, label="insulation")
		max=np.max(np.array(ins_mid['insulationScore']).ravel())
		min=np.min(np.array(ins_mid['insulationScore']))
		ax.set_xlim(self.real_start,self.real_end)
		ax.set_xticks([])
		ax.set_yticks([min,max])
		ax.set_yticklabels([min,max],fontsize=3)
		for loc in ['left','right','bottom']:
			ax.spines[loc].set_linewidth(0.5)
			ax.spines[loc].set_color('k')
		ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
		ax.text(-0.01,0.5,'insulations',horizontalalignment='right',verticalalignment='center',transform=ax.transAxes,fontsize=4)
		return line

	def cut_insulation(self):
		ins_df=pd.read_table(self.insulation,header=0,index_col=None,encoding="utf-8")
		ins_df=ins_df.fillna(0)
		ins_df=ins_df.loc[:,['start','end','insulationScore']]
		ins_df['mid']=(ins_df['start']+ins_df['end'])/2
		ins_df=ins_df.loc[:,["mid","insulationScore"]]
		ins_mid=ins_df[(ins_df['mid']>(self.real_start-1))&(ins_df['mid']<(self.real_end+1))]
		ins_mid=ins_mid.reset_index(drop=True)
		return ins_mid

	def prepare_TAD_lay_out(self):
		'''
			prepare the plot layout
			what to consider is how to set fig width 
			and the ratio between heatmap and boundary and peak
			contact heatmap : 5 grid
			gene Denesity map : 1 grid
			boundaries map : 3 grid 
			insulation map : 3 grid
		'''
		grids=2
		grid_layouts=[2,]
		if self.show_map:
			grids+=2;grid_layouts.append(grids)
		if self.show_boundaries: 
			grids+=1;grid_layouts.append(grids)
		if self.show_insulation:
			grids+=1;grid_layouts.append(grids)
		self.grids=grids
		self.tad_grid_layout=grid_layouts[::-1]
		G=gridspec.GridSpec(grids,1)
		G.update(hspace=0)
		self.Gs=G

	def draw_TAD(self,fig):
		G=self.Gs
		l=self.tad_grid_layout
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		print('TAD layout......')
		print(l)
		grids=self.grids
		print('grids.....')
		print(grids)
		plt.suptitle(self.title,fontsize=12)
		#----------------draw trimap ----------------------#
		if self.show_map:
			#ax1 = plt.subplot(G[l.pop():l[-1],:],axisbg='white')
			gridsSub = gridspec.GridSpecFromSubplotSpec(1,1, subplot_spec=G[l.pop():l[-1],:], wspace=0, hspace=0)
			ax1=fig.add_subplot(gridsSub[:,:],axisbg='white')
			self.sample_mat.tri_map(ax1,self.sample_mat.all_mat,False)
		#----------------draw boundary --------------------#
		if self.show_boundaries:
			#ax2=plt.subplot(G[l.pop():l[-1],:],axisbg='white')
			gridsSub2 = gridspec.GridSpecFromSubplotSpec(1,1, subplot_spec=G[l.pop():l[-1],:], wspace=0, hspace=0)
			ax2=fig.add_subplot(gridsSub2[:,:],axisbg='white')
			self.draw_bounds(ax2)
		#---------------draw insulation -------------------#
		if self.show_insulation:
			#ax3=plt.subplot(G[l.pop():l[-1],:],axisbg='white')
			gridsSub3 = gridspec.GridSpecFromSubplotSpec(1,1, subplot_spec=G[l.pop():l[-1],:], wspace=0, hspace=0)
			ax3=fig.add_subplot(gridsSub3[:,:],axisbg='white')
			self.draw_insulation(ax3)

def adjust_fig_size(show_map=True,show_insulation=True,show_boundaries=True,start=0,end=0):
	grids=2
	grid_layouts=[2,]
	if show_map:
		grids+=3
	if show_boundaries:
		grids+=1
	if show_insulation:
		grids+=2
	height_unit=0.25
	width_unit=0.25
	fig_height=grids*height_unit
	fig_width=width_unit*((end-start)/1000000)
	return fig_width,fig_height

class TestdekkerTAD(unittest.TestCase):
	def test_sample(self):
		dekkerTad=TAD_insulation('DM','chr19',1000000,61342430)
		dekkerTad.mat_path='{}/../Example/single_insulation_TAD/data/DM1.mm9_chr19.mm9_chr19.matrix'.format(Bin)
		dekkerTad.boundaries='{}/../Example/single_insulation_TAD/data/DM1_chr19_immeran.boundaries'.format(Bin)
		dekkerTad.insulation='{}/../Example/single_insulation_TAD/data/DM1_chr19_immeran.insulation'.format(Bin)
		dekkerTad.resolution=40000
		dekkerTad.extend=''
		dekkerTad.region=''
		dekkerTad.outdir="{}/../Example/single_insulation_TAD/".format(Bin)
		fig = plt.figure(figsize=(16,2),dpi=1200)
		dekkerTad.prepare_TAD_lay_out()
		dekkerTad.draw_TAD(fig)
		out_fig="{}/{}_{}_{}_{}_boundaries.pdf".format(dekkerTad.outdir,dekkerTad.name,dekkerTad.chr,dekkerTad.start,dekkerTad.end)
		fig.savefig(out_fig)
		print(out_fig)

if __name__ == '__main__':
	unittest.main()
