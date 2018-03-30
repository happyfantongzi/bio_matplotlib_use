'''
	this is used cut matrix and draw map
	1.cut matrix
	2. tranform matrix
	need:
		matrix  or  matrix file
		resolution
		genomic region
	method:
		cut_mat
		normlize_mat
		heatmap
		Triangular
		hexagonal
'''
import math,re,sys,os
import numpy as np
old_settings = np.seterr(all='ignore')
import pandas as pd
pd.set_option('display.precision', 3)
import unittest
Bin = os.path.abspath(os.path.dirname(__file__))
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['xtick.direction'] = 'out' 
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['font.size']=8
from matplotlib import pyplot as plt
plt.style.use('ggplot')

__author__='fanxuning'
__mail__='xuningfan@genome.cn'

class subMatrix():
	def __init__(self,chr,start,end,resolution,focusRegions):
		self.sample=''
		self.mat_path=''
		self.mat=pd.DataFrame()
		self.need_to_be_cut=True
		self.chr=chr
		self.start=start
		self.end=end
		self.resolution=resolution
		self.focusRegion=focusRegions
		self.showFocusRegion=False
		self.showScale=True
		#self.showRatio=0.1
		self.showRatio=0.04
		self.up_percent=25
		self.down_percent=95
		self.colormap=plt.get_cmap('YlOrBr')
		self.log_transform=False
		self.need_to_normalized=False
		self.norm_method='z_score_whole_genome'
		self.vmin=0
		self.vmax=0
		self.up_down_together=False
		self.extract_percent=90
		self.up=0;self.down=0

	@property
	def reStart(self):
		start=int(self.start/self.resolution)
		return start

	@property
	def genome_start(self): 
		genome_start=self.reStart*self.resolution
		return genome_start

	@property
	def reEnd(self):
		end=math.ceil(self.end/self.resolution)
		return end

	@property
	def genome_end(self):
		genome_end=self.reEnd*self.resolution
		return genome_end

	@property
	def all_mat(self):
		if not self.mat_path=='':
			mat=pd.read_table(self.mat_path,header=0,index_col=0)
			return mat
		elif not self.mat.shape[0]==0:
			mat=self.mat
			return mat
		else :
			sys.exit('please input matrix or matrix file!!!')

	def cut_mat(self,mat):
		#self.all_mat
		sub_mat=mat.iloc[self.reStart:self.reEnd,self.reStart:self.reEnd]
		sub_mat=sub_mat.reset_index(drop=True)
		return sub_mat

	def normalize_matrix(self,raw_mat,type='z_score_whole_genome'):
		all_norm_mat=pd.DataFrame()
		desc = pd.DataFrame(np.ravel(raw_mat)).describe(include=[np.number])
		if type=='z_score_whole_genome':
			all_norm_mat = (raw_mat - desc.ix['mean',0])/desc.ix['std',0]
		if type=='min_max':
			all_norm_mat=(raw_mat-desc.ix['min',0])/(desc.ix['max',0]-desc.ix['min',0])
		if type=='z_score_by_distance':
			curmat=np.array(raw_mat)
			mat_len=curmat.shape[0]
			expect=np.zeros([mat_len,mat_len])
			sd=np.zeros([mat_len,mat_len])
			for i in range(mat_len):
				curlist=curmat[range(0,mat_len-i),range(i,mat_len)]
				curmean=np.mean(curlist)
				curstd=np.std(curlist)
				expect[range(0,mat_len-i),range(i,mat_len)]=curmean
				sd[range(0,mat_len-i),range(i,mat_len)]=curstd
			z_score_mat=(curmat-expect)/sd
			z_score_mat[np.isnan(z_score_mat)]=0
			z_score_mat[np.isinf(z_score_mat)]=0
			all_norm_mat=z_score_mat
		return all_norm_mat

	def Triangular(self,raw_mat):
		mat=pd.DataFrame()
		if self.need_to_normalized:
			all_norm_mat=self.normalize_matrix(raw_mat,type=self.norm_method)
			mat=self.cut_mat(all_norm_mat) if self.need_to_be_cut else all_norm_cut
		else:
			mat=self.cut_mat(raw_mat) if self.need_to_be_cut else raw_mat 
		mat=np.array(mat)
		n=0;length=len(mat)
		tri_mat=np.zeros([length,length*2])
		tri_mat[tri_mat==0]=np.nan
		for i in range(length):
			curl=np.array(np.diag(mat,i))
			tri_mat[i,n:(length*2-n)]=curl.repeat(2)
			n+=1
		return tri_mat 

	def hexagonal(self,raw_mat):
		mat=pd.DataFrame()
		if self.need_to_normalized:
			all_norm_mat=self.normalize_matrix(type=self.norm_method)
			mat=self.cut_mat(all_norm_mat)
		else:
			mat=self.cut_mat(raw_mat)
		mat=np.array(mat)
		length=len(mat)
		hex_mat=np.zeros([length*2,length*2])
		n=length-1
		for i in range(-length+1,length):
			curl=np.array(np.diag(mat,i)).repeat(2)
			shape=curl.shape[0]
			if i<=0:
				hex_mat[i+length-1,n:n+shape]=curl
				hex_mat[length-i-1,n:n+shape]=curl
				n=n-1
		return hex_mat

	def adjust_mat(self,raw_mat,type,zoom_out=False): 
		'''
			two ways:
				the first way : shortend the difference,to show TAD
				the second way : zoom up the difference of two TAD
		'''
		up_percent=self.up_percent
		down_percent=self.down_percent
		triMat=self.Triangular(raw_mat)
		triArr=np.array(triMat)
		if self.log_transform:
			triArr=np.log(triMat)
		if self.up_down_together:
			tmp_trimat=np.abs(triArr)
			xmax=np.nanpercentile(tmp_trimat,self.extract_percent)
			self.vamx=xmax;self.vmin=-self.vmax
		else:
			xmax=np.nanpercentile(triArr,down_percent)
			xmin=np.nanpercentile(triArr,up_percent)
			if type=='original':
				self.vmax=xmax;self.vmin=xmin
			else:
				self.vmax=max(abs(xmax),abs(xmin));self.vmin=-self.vmax
		if zoom_out:
			triArr[(triArr>xmin)&(triArr<xmax)]=0
		return triArr

	def colorbar(self,ax,im,type):
		print('{} trimap colorbar'.format(self.sample))
		print('vmin:{}\tvmax:{}'.format(self.vmin,self.vmax))
		print('vmin:{}\tvmax:{}'.format(round(self.vmin),round(self.vmax)))
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		if type=='trimap':
			axins1 = inset_axes(ax, width="1%", height="50%",loc=3, bbox_to_anchor=(0, 0.2, 0.5, 1), bbox_transform=ax.transAxes,borderpad=0)
			cbar=plt.colorbar(im, cax=axins1, orientation='vertical',ticks=[round(self.vmin),int(self.vmax)])
		else:
			axins1 = inset_axes(ax, width="100%", height="100%",loc=6, bbox_to_anchor=(1.03, 0.4, 0.04, 0.36), bbox_transform=ax.transAxes,borderpad=0)
			cbar=plt.colorbar(im, cax=axins1, orientation='vertical',ticks=[math.ceil(self.vmin),math.floor(self.vmax)])
		return cbar

	def add_colorbar(self,im, aspect=20, pad_fraction=0.5, **kwargs):
		from mpl_toolkits import axes_grid1
		from mpl_toolkits.axes_grid1.inset_locator import inset_axes
		divider= axes_grid1.make_axes_locatable(im.axes)
		width= axes_grid1.axes_size.AxesY(im.axes, aspect=1/aspect)
		pad=axes_grid1.axes_size.Fraction(pad_fraction, width)
		cax= divider.append_axes("right", size=width, pad=pad)
		plt.sca(current_ax)
		return im.axes.figure.colorbar(im, cax=cax, **kwargs)


	def plotHeatmap(self,raw_mat,ax):
		'''
			1.remove nan 
			2. colorbar set....
			still to be done.....
		'''
		for loc in ['left','top','right','bottom']:
			ax.spines[loc].set_linewidth(0.5)
			ax.spines[loc].set_color('k')
		ax.tick_params(bottom ='on',top='off',left='on',right='off')
		up =np.nanpercentile(raw_mat,self.up_percent)
		down= np.nanpercentile(raw_mat,self.down_percent)
		self.up=up;self.down=down
		norm = mpl.colors.Normalize(up, down)
		cax =ax.pcolormesh(np.array(raw_mat),cmap=self.colormap,norm=norm)
		#cax = ax.matshow(np.array(raw_mat),cmap=self.colormap,clim=(up,down))
		bins=self.reEnd-self.reStart
		ticks=list(np.linspace(0,bins,5))
		ax.set_xticks(ticks);ax.set_yticks(ticks)
		#ticks_labs=list(np.linspace(self.genome_start,self.genome_end,5))
		ticks_labs=list(np.linspace(self.genome_start/1000000,self.genome_end/1000000,5))
		#tick_labels=['{:,}'.format(int(x)) for x in ticks_labs]
		tick_labels=[int(x) for x in ticks_labs]
		#ax.set_xticklabels(tick_labels,fontsize=3,rotation=45);ax.set_yticklabels(tick_labels,fontsize=3,rotation=45)
		ax.set_xticklabels(tick_labels);ax.set_yticklabels(tick_labels)
		bbox_to_anchor=(1.03, 0.4, 0.04, 0.36)
		loc=6;self.colorbar(ax,cax,'heatmap')
		if not self.sample=='':
			ax.set_title('{}'.format(self.sample),fontsize=15)
		ax.set_xlim(0,bins)
		ax.set_ylim(0,bins)
		ax.grid(False)

	def mark_out_interested_region(self,ax,region_list):
		'''
		two heatmap type :
			1.tri map make out :
					just show genome region , annotating......
			2.heatmap make out :
					a. genome region ,it will like a strip...
					b.interactive  bin pairs like a circle
		mutiple makers to be choose:
			1.arrow
			2.circle
			3.rectangle etc.
		steps and notice:
			1.set(x,y)
			2.compare with bin involved
			3.show interface.....
 		a lot can  be done.....
		'''
		from matplotlib.patches import Circle
		for region in region_list:
			region = Circle((x,y),cisize, fc='none',ec='black')
			ax.add_patch(region)


	def tri_map(self,ax,raw_mat,flag):
		ax.tick_params(top='off',bottom='on',left='off',right='off')
		trimat=self.adjust_mat(raw_mat,'original',zoom_out=False)
		for loc in ['left','right','top']:
			ax.spines[loc].set_visible(False)
		#-------draw tri map -------------------------------------#
		with np.errstate(divide='ignore'):
			#im=ax.imshow(trimat,origin="lower",cmap=self.colormap,interpolation="nearest",aspect='auto')
			im=ax.matshow(trimat,cmap=self.colormap,clim=(self.vmin,self.vmax),interpolation="nearest",aspect='auto')
			self.colormap.set_bad('white')
		#-------set ticks------------------------------------------#
		ax.set(yticks=[])
		ax.tick_params(direction='out',pad=5)
		length=len(trimat)*2
		ticks=np.linspace(0,length,5)
		ticklabs=np.linspace(self.reStart*self.resolution,self.reEnd*self.resolution,5)
		tkls=['{:,}'.format(int(i)) for i in ticklabs ]
		#-------set ylim-------------------------------------------#
		ax.set_ylim(0,self.showRatio*length)
		if flag:
			ax.xaxis.tick_bottom()
			ax.spines['bottom'].set_position(('outward',2))
			ax.spines['bottom'].set_linewidth(0.5)
			ax.spines['bottom'].set_color('k')
			print('ticks....')
			print(ticks)
			print('ticklables...')
			print(tkls)
			ax.set_xticks(ticks)
			ax.set_xticklabels(tkls,fontsize=3)
		else:
			ax.set_xticks([])
		#--------focus genome region eclipse -------------------------------#
		self.colorbar(ax,im,'trimap')
		if self.showFocusRegion:
			print('draw focus region')
		ax.grid(False)
		return im
	

class TestMatrix(unittest.TestCase):
	def test_cut_matrix(self):
		data_dir='{}/../Example/ComputMatrix/data/'.format(Bin)
		path='{}/DM.mm9_chr19.mm9_chr19.matrix'.format(data_dir)
		chr='chr19'
		resolution=40000;focusRegion=0
		mat=subMatrix(chr,0,61342430,resolution,focusRegion)
		mat.sample='DM'
		mat.colormap=plt.get_cmap('YlOrRd')
		mat.mat_path=path
		fig=plt.figure(figsize=(16,10))
		ax=fig.add_axes([0.05,0.85,0.8,0.05],axisbg='w')
		mat.tri_map(ax,mat.all_mat,True)
		outdir='{}/../Example/ComputMatrix/'.format(Bin)
		print('out  trimap....')
		print("{}/unittest.chr19_tri_map.pdf".format(outdir))
		fig.savefig('{}/unittest.chr19_tri_map.pdf'.format(outdir))
		fig2=plt.figure(figsize=(4,4))
		ax=fig2.add_axes([0.1,0.1,0.8,0.8],axisbg='w')
		raw_mat=mat.cut_mat(mat.all_mat)
		mat.plotHeatmap(raw_mat,ax)
		print('out heatmap.....')
		print('{}/chr19_heatmap.pdf'.format(outdir))
		fig2.savefig('{}/chr19_heatmap.pdf'.format(outdir))

if __name__ == '__main__':
	unittest.main()

