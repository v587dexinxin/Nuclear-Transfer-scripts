# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:48:23 2019

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
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['skyblue' , 'k' , 'yellow'])
#my_cmap.set_bad('#2672a1')

my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'violet'])
my_cmap.set_bad('#2672a1')

def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'transcript_id',  'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'transcript':
            gene_name = i.strip().split('gene_name')[1].split('\"')[1]
            gene_id = i.strip().split('gene_id')[1].split('\"')[1]
            transcript_id = i.strip().split('transcript_id')[1].split('\"')[1]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , transcript_id , gene_name , chro , strand , start , end))
    gtf.close()
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf
    
def plot_heatmap(matrix):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    vmax = np.percentile(matrix,95)
    vmin = np.percentile(matrix,5)
    im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    plt.colorbar(im)
    x = ['CCS','NT5','NT6','fESC']
    numbers = {'1':795 , '2':3155}

    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 30)
    ax.set_ylabel('\n'.join(['cluster' + str(int(i)) + ':' + str(numbers[str(int(i))]) for i in [1,2]]) , fontsize = 20)
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('Zscore ,gene numbers:' + str(len(matrix)),fontsize = 30)
    return fig

def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()

geneFil = 'H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt'
gtfFil = 'G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf'

gtf_type = np.dtype({'names':['gene_id' , 'transcript_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S64' , 'S8' , 'S8' , np.int , np.int]})

gene_type = np.dtype({'names':['gene_name' , 'chr' , 'starnd' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                      'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

id_type = np.dtype({'names':['transcript_id'],
                            'formats':['S64']})

gtf = Load_gtf(gtfFil)
geneData = np.loadtxt(geneFil , dtype = gene_type , skiprows = 1)



transcript_id = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\promoter\\union_cluster2_promoter.txt' , 
                           dtype= id_type, usecols = (12) , skiprows=1)

genes =[]

n = 0
for i in transcript_id:
    mask = (gtf['transcript_id'] == i['transcript_id'])
    overlap = gtf[mask]
    if overlap.size != 0:
        gene_name = overlap['gene_name']
        gene = geneData[geneData['gene_name'] == gene_name][0]
        genes.append(gene)
    else:
        n += 1
        
genes = np.array(genes , dtype = gene_type)

CCS_gene = []
for i in genes:
    if (i['CCS'] / i['NT5'] > 1.5) and (i['CCS'] / i['NT6'] > 1.5) and (i['CCS'] / i['fESC'] > 1.5):
        CCS_gene.append(i)
CCS_gene = np.array(CCS_gene , dtype = gene_type)
        
ESC_gene = []
for i in genes:
    if (i['NT5'] / i['CCS'] > 1.5) and (i['NT6'] / i['CCS'] > 1.5) and (i['fESC'] / i['CCS'] > 1.5):
        ESC_gene.append(i)
ESC_gene = np.array(ESC_gene , dtype = gene_type)

w = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\union_cluster1_pro_geneID.txt' , 'w')
for i in genes:
    gene_name = i['gene_name']
    g_id = gtf[gtf['gene_name'] == gene_name][0]['gene_id']
    w.writelines(g_id + '\n')
w.close()

w = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\union_cluster3_pro_gene.txt' , 'w')
for i in ESC_gene:
    gene_name = i['gene_name']
    g_id = gtf[gtf['gene_name'] == gene_name][0]['gene_name']
    w.writelines(g_id + '\n')
w.close()   

##Get_FPKM_Matrix
matrix = np.zeros((len(CCS_gene) + len(ESC_gene) , 4))

for i in range(len(CCS_gene)):
    matrix[i,0] = np.log2(CCS_gene[i]['CCS'] + 1)
    matrix[i,1] = np.log2(CCS_gene[i]['NT5'] + 1)
    matrix[i,2] = np.log2(CCS_gene[i]['NT6'] + 1)
    matrix[i,3] = np.log2(CCS_gene[i]['fESC'] + 1)
    
for i in range(len(CCS_gene) , len(CCS_gene) + len(ESC_gene)):
    matrix[i,0] = np.log2(ESC_gene[i - len(CCS_gene)]['CCS'] + 1)
    matrix[i,1] = np.log2(ESC_gene[i - len(CCS_gene)]['NT5'] + 1)
    matrix[i,2] = np.log2(ESC_gene[i - len(CCS_gene)]['NT6'] + 1)
    matrix[i,3] = np.log2(ESC_gene[i - len(CCS_gene)]['fESC'] + 1)

matrix_0 = Z_score(matrix)

#Plot    
fig = plot_heatmap(matrix_0)    
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\Selected_gene_heatmap_all_Zscore_bwr.pdf')
    
    