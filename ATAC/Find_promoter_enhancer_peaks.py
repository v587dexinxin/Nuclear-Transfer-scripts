# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 08:19:32 2019

@author: han-luo
"""

from __future__ import division
import math
import numpy as np
import csv
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

def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})
    
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene


def get_promoter_pos(union_gene):
    pro_type = np.dtype({'names':['gene_name' , 'chr' , 'start' , 'end' , 'strand'],
                      'formats':['<S64' , '<S8' , np.int , np.int , '<S2']})
    promoter = []
    for gene in union_gene:
        if gene['strand'] == '+':
            promoter.append((gene['gene_name'] , gene['chr'] , gene['start'] - 1000 , gene['start'] + 1000 , gene['strand']))
            promoter.append((gene['gene_name'] , gene['chr'] , gene['end'] - 1000 , gene['end'] + 1000 , gene['strand']))
        else:
            promoter.append((gene['gene_name'] , gene['chr'] , gene['start'] - 1000 , gene['start'] + 1000 , gene['strand']))
            promoter.append((gene['gene_name'] , gene['chr'] , gene['end'] - 1000 , gene['end'] + 1000 , gene['strand']))
    promoter = np.array(promoter , dtype = pro_type)
    return promoter
    
def annotation_peak(PeakFil , proFil, enFil):
    '''
    '''
    
    peaktype = np.dtype({'names':['chr','start','end'],
                      'formats':['<S8' , np.int , np.int]})
    peaks = np.loadtxt(PeakFil , usecols = (0 , 1 , 2) , dtype = peaktype)
    
    w1 = open(proFil,'w')
    w2 = open(enFil,'w')
    n = 0
    for g in set(peaks['chr']):
        tmp_1 = peaks[peaks['chr'] == g]
        tmp_2 = promoter[promoter['chr'] == g]
        for p in tmp_1:
            mask = (tmp_2['start'] <= p['end']) & (tmp_2['end'] > p['start'])
            overlap = tmp_2[mask]
            p = np.array(list(p),dtype = np.str)
            if overlap.size != 0:
                n += 1
                if overlap.size > 1:
                    print 'peaks:' , p , overlap
                w1.writelines('\t'.join(p) + '\t' + ','.join([i['gene_name'] for i in overlap]) + '\n')
            else:
                w2.writelines('\t'.join(p) + '\n')
    
    w1.close()
    w2.close()
    return n 

#------------------------get_union_gene---------------------------------
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')  

#------------------------get_promoter_pos-------------------------------
promoter = get_promoter_pos(union_gene)

PeakFil = 'H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union\\union_cluster3.bed'
proFil = 'H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union\\union_cluster3_Promoter.bed'
enFil = 'H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union\\union_cluster3_Enhancer.bed'
n = annotation_peak(PeakFil , proFil , enFil)



