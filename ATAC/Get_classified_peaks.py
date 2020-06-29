# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:04:16 2019

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



def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a
    
union_type = np.dtype({'names':['chr' , 'start' , 'end' , 'groups'],
                       'formats':['S8' , np.int , np.int , 'S64']})

union = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\union_Kmeans4sortedRegions.bed' , dtype = union_type , skiprows = 1 , usecols = (0,1,2,-1))

cluster1 = [] ; cluster2 = [] ; cluster3 = []

for i in union:
    if i['groups'] == 'cluster_1':
        cluster1.append(i)
    elif i['groups'] == 'cluster_4':
        cluster2.append(i)
    else:
        cluster3.append(i)

cluster1 = Sort(cluster1 , 0 , 1)
cluster2 = Sort(cluster2 , 0 , 1)        
cluster3 = Sort(cluster3 , 0 , 1) 

w = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\union_cluster2.bed' , 'w')
for i in cluster2:
    w.writelines('\t'.join(['chr' + i[0] , str(i[1]) , str(i[2])]) + '\n')
w.close()

                 