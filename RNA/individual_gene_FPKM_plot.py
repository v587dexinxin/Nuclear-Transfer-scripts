# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 17:51:25 2019

@author: han-luo
"""
from __future__ import division
import math
import numpy as np
import csv , copy
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Our Own Color Map
#my_cmap = plt.get_cmap('bwr')
#my_cmap.set_bad('#2672a1')

my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'yellow'])
my_cmap.set_bad('#2672a1')

gene_type = ({'names':['gene_name'],
               'formats':['S64']})
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()

Runx2 = np.log2(np.array([14.165874 , 20.797432 , 21.22489 , 0.193843,	0.110375,	0.282659	,0.487601,	0.256296,	1.792519	,2.554883,	0.795944,	0.58317	,0.472683]) + 1)
Oct4 = np.log2(np.array([0.334727	,0.460841	,0.263261	,642.058838,	609.059265	,639.192383	,605.897644	,514.524292	,527.433472,	460.26178	,589.205139	,536.029175,	563.454651]) + 1)
Nanog = np.log2(np.array([0.630105	,0.403784	,0.432458,	518.177673	,571.663574	,568.905579	,434.392242,	480.410919,	313.998688	,280.508057,	421.953156,	527.921021	,517.037842]) + 1)
Sox2 = np.log2(np.array([0.155137 , 0.188073	, 0.159787 , 205.959229 , 150.405426 , 129.225784 , 174.896729	 , 194.115784	 , 214.367081 , 158.780472 , 171.021973 , 213.716202]) + 1)
Parm1 = np.log2(np.array([467.20163 , 289.621552 , 331.0271 , 1.63998 , 0.773053 , 0.653098 , 0.942214 , 1.424455 , 1.701045 , 1.538032 , 1.626868 , 1.370357]) + 1)


Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
fig = plt.figure(figsize = (8, 8))
ax = fig.add_axes([Left  , HB , width , HH])
ax.scatter([1,1,1] , Oct4[:3] ,s = 700)
ax.scatter([2,2,2] , Oct4[6:9] ,s = 700)
ax.scatter([3,3,3] , Oct4[10:],s = 700)
ax.scatter([4,4,4] , Oct4[3:6],s = 700)
ax.scatter([1,2,3,4] , [Oct4[:3].mean() , Oct4[6:9].mean() , Oct4[10:].mean() , Oct4[3:6].mean()] , marker = '_' , s = 1000 , c = 'black')
ax.set_xlim((0.5,5))
ax.set_ylim((0,10))
ax.set_xticks([1,2,3,4])
ax.set_xticklabels(['CCS' , 'NT5' , 'NT6' , 'fESC'] ,fontsize = 20)
ax.set_ylabel('np.log2(FPKM+1)' , fontsize = 30)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\Nanog_FPKM.pdf')

Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
fig = plt.figure(figsize = (8, 8))
ax = fig.add_axes([Left  , HB , width , HH])
ax.scatter([1,1,1] , Sox2[:3] ,s = 700)
ax.scatter([2,2,2] , Sox2[3:6] ,s = 700)
ax.scatter([3,3,3] , Sox2[6:9],s = 700)
ax.scatter([4,4,4] , Sox2[9:12],s = 700)
ax.scatter([1,2,3,4] , [Sox2[:3].mean() , Sox2[6:9].mean() , Sox2[10:].mean() , Sox2[3:6].mean()] , marker = '_' , s = 1000 , c = 'black')
ax.set_xlim((0.5,5))
ax.set_ylim((0,10))
ax.set_xticks([1,2,3,4])
ax.set_xticklabels(['CCS' , 'NT5' , 'NT6' , 'fESC'] ,fontsize = 20)
ax.set_ylabel('np.log2(FPKM+1)' , fontsize = 30)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\H_Sox2_FPKM.pdf')