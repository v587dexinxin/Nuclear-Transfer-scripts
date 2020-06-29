# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 15:20:20 2020

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
               
def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
            union.append(i)
    union = np.array(union , dtype = union_type)
    return union
    
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()


union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')


Magea2 = np.log2(np.array(list(union_gene[union_gene['gene_name'] == 'Magea2'][0])[5:]) + 1)
Magea3 = np.log2(np.array(list(union_gene[union_gene['gene_name'] == 'Magea3'][0])[5:]) + 1)
Magea5 = np.log2(np.array(list(union_gene[union_gene['gene_name'] == 'Magea5'][0])[5:]) + 1)
Magea6 = np.log2(np.array(list(union_gene[union_gene['gene_name'] == 'Magea6'][0])[5:]) + 1)
Magea8 = np.log2(np.array(list(union_gene[union_gene['gene_name'] == 'Magea8'][0])[5:]) + 1)

Magea8 = np.log2(np.array(list(union_gene[union_gene['gene_name'] == 'Magea8'][0])[5:]) + 1)
Magea8 = np.log2(np.array(list(union_gene[union_gene['gene_name'] == 'Magea8'][0])[5:]) + 1)





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