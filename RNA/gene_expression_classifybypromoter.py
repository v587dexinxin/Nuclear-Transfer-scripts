# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 08:56:27 2019

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
    return union_gene

def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix


def get_indexed_matrix(gene_name , zscore):
    matrix = np.zeros((len(gene_name) , 4))
    for i in range(len(gene_name)):
        gene = union_gene[union_gene['gene_name'] == gene_name[i]][0]
        matrix[i][0] = np.log2(gene['CCS'] + 1)
        matrix[i][1] = np.log2(gene['NT5'] + 1)
        matrix[i][2] = np.log2(gene['NT6'] + 1)
        matrix[i][3] = np.log2(gene['fESC'] + 1)
    if zscore == True:
        matrix = Z_score(matrix)
    else:
        pass
    index = np.arange(len(gene_name))
    index = index.reshape((len(gene_name) , 1))
    matrix = np.hstack((matrix , index))
    return matrix

def plot_heatmap(matrix):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    vmax = np.percentile(matrix[:,:-1],95)
    vmin = np.percentile(matrix[:,:-1],5)
    im = ax.imshow(matrix[:,:-1],vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    plt.colorbar(im)
    x = ['CCS','NT5','NT6','fESC']
    numbers = {'0':131 , '1':296}

    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 30)
    ax.set_ylabel('\n'.join(['cluster' + str(int(i)) + ':' + str(numbers[str(int(i))]) for i in [0,1]]) , fontsize = 20)
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('log2(FPKM + 1) ,gene numbers:' + str(len(matrix)),fontsize = 30)
    return fig

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
 
 
               
##-----------------------------get_diff_gene_names--------------------------------------------
diff_gene = np.loadtxt('H:\\Workspace_New\\data\\RNA\\diff_expression\\diff_gene_matrix_np.log2(FPKM+1)_classify2.txt' , dtype = gene_type , usecols = (0) , skiprows = 1)
diff_gene = diff_gene['gene_name']  

#------------------------------get_union_gene----------------------------------------------
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')        

##-----------------------------get_gene_name--------------------------------------------------   

CCS_name = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\classify\\CCS_noNT_nofESC_Promoter.bed' , dtype = gene_type , usecols = (-1))
CCS_name = CCS_name['gene_name']
ES_name = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\classify\\noCCS_NT_fESC_Promoter.bed' , dtype = gene_type , usecols = (-1))
ES_name = ES_name['gene_name']

gene_name = []
for i in CCS_name:
    if i not in gene_name:
        gene_name.append(i)

for i in ES_name:
    if i not in gene_name:
        gene_name.append(i)

gene_names = []
for i in gene_name:
    i = set(i.split(','))
    for j in i:
        if j in diff_gene:
           gene_names.append(j)
        else:
            pass

len(gene_names)


#-----------------------------get_plot_matrix-------------------------------------------
matrix_0 = get_indexed_matrix(gene_names , False)

#--------------------------------plot-----------------------------------------
fig = plot_heatmap(matrix_0)
outFil = 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\Diff_promoter_gene_expression_heatmap_np.log2(FPKM+1).pdf'    


run_Plot(fig , outFil)
        