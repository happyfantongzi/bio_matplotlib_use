import argparse
import re,sys,os
import numpy as np
np.seterr(invalid='ignore')
import pandas as pd
import unittest

__author__='fanxuning'
__mail__='xuningfan@genome.cn'

class Bounds_classifier: 
	def __init__(self,sam1,sam2):
		self.sam1=sam1
		self.sam2=sam2
		self.bound1=pd.DataFrame()
		self.bound2=pd.DataFrame()
		self.bound1_path=''
		self.bound2_path=''
		self.interval=80000
		self.resolution=40000

	def extract_common_border(self):
		'''
			extract_common_border.....
		'''
		control=self.bound1['start']
		treat=self.bound2['start']
		common_list=list(set(control) & set(treat))
		qualified_pairs=dict(zip(common_list,common_list))
		con_uniqs=list(set(control).difference(set(treat)))
		treat_uniqs=list(set(treat).difference(set(control)))
		from itertools import product
		reflect_pair=list(product(con_uniqs,treat_uniqs))
		for item in reflect_pair:
			inter=abs(item[0]-item[1])
			if inter <= self.interval:
				qualified_pairs[item[0]]=item[1]
		return qualified_pairs

	def classify_boundaries(self):
		qualified_pairs=self.extract_common_border()
		self.bound1['type']=self.bound1['start'].apply(
				lambda x: 'overlaping'
				if x in qualified_pairs.keys()
				else 'dispearing'
				)
		self.bound1['start']=self.bound1['start'].apply(
				lambda x:(x+qualified_pairs[x])/2 if x in qualified_pairs.keys()
				else x
				)
		self.bound2['type']=self.bound2['start'].apply(
				lambda x: 'overlaping'
				if x in qualified_pairs.values() 
				else 'reforming'
				)
		reforming=self.bound2[self.bound2['type']=='reforming']
		reforming=reforming.reset_index(drop=True)
		classify_bounds=pd.concat([self.bound1,reforming])
		classify_bounds=classify_bounds.sort_values(by='start')
		classify_bounds=classify_bounds.reset_index(drop=True)
		return classify_bounds

class Test_bound_classifier(unittest.TestCase): 
	def test_classifier(self):
		print('test  boundaries classifier......')

if __name__=='__main__':
	unittest.main()
