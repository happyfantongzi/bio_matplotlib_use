#! /usr/bin/env python3
import re,sys,os
import argparse

'''
	this is just for temp use 
	In fact:
		what's i wanna do .....
			is build a systermic Hic advance analysis
			so, i will keep up to complete it and update it
			After......
			i will add machine learning model ......
'''
__author__='fanxuning'
__mail__='xuningfan@genome.cn'

def read_Delta_hic_options():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}\n'.format(__author__,__mail__))
	parser.add_argument('-gr','--genomeSite',help='geneome site',dest='genome_region',type=str,required=True)
	'''
			sam1	sam2
		mat	mat1	mat2
		bound	bound1	bound2
		listbam	listbam1	listbam2
		gene	anno	anno
	'''


	args=parser.parse_args()


def main():


if __name__=='__main__':
	main()
