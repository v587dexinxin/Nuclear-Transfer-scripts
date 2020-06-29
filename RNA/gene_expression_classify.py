# -*- coding: utf-8 -*-
"""
Created on Mon Oct 07 14:45:12 2019

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

def get_gene_exp(union_Fil , genename):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(union_Fil , dtype = union_type , skiprows = 1)
    gene = []
    for i in genename:
        g = union_gene[union_gene['gene_name'] == i]
        gene.append(g[0])
    gene = np.array(gene , dtype = union_type)
    return gene

def get_union_gene_name(cell):
    gene_name = []
    for c in cell:
        c1 = c.split("_")[0] ; c2 = c.split("_")[1]
        gene_type = ({'names':['gene_id' , 'baseMean' , 'log2FoldChange' , 'lfcSE' , 'stat' , 'pvalue' , 'padj' , 'gene_name' , 'chr' , 'strand' , 
                          'start' , 'end' , c1  , c2 ],
                     'formats':['S64' , np.float , np.float , np.float , np.float , np.float , np.float , 'S64' , 'S8' , 'S8' , np.int , np.int , 
                                np.float , np.float]})
        diffData = np.loadtxt('D:\\Workspace_New\\data\\RNA\\diff_expression\\basemean_500_fc_1.5\\Filtered_' + c1 + '_' + c2 + '.txt' , skiprows = 1 , 
                             dtype = gene_type)
        for i in diffData:
            gene_name.append(i['gene_name'])
        
    gene_name = list(set(gene_name))
    return gene_name

def get_common2_gene_name(Fil1 , Fil2):
    diffData1 = np.loadtxt(Fil1 , skiprows = 1 , usecols = (7) , dtype = 'S64')
    diffData2 = np.loadtxt(Fil2 , skiprows = 1 , usecols = (7) , dtype = 'S64')
    intersection = set(diffData1).intersection(set(diffData2))
    
    return intersection


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

def K_means_cluster(matrix_0,matrix,n , reverse):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1] , reverse=reverse))
    return matrix_1

    
def plot_heatmap(matrix):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    vmax = np.percentile(matrix[:,:-2],95)
    vmin = np.percentile(matrix[:,:-2],5)
    im = ax.imshow(matrix[:,:-2],vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    plt.colorbar(im)
    x = ['CCS','NT5','NT6','fESC']
    numbers = {}
    for i in set(matrix[:,-1]):
        numbers[str(int(i))] = 0
        for j in matrix:
            if i == j[-1]:
                numbers[str(int(i))] += 1
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 30)
    ax.set_ylabel('\n'.join(['cluster' + str(int(i)) + ':' + str(numbers[str(int(i))]) for i in [1,0]]) , fontsize = 20)
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('log2(FPKM + 1) ,gene numbers:' + str(len(matrix)),fontsize = 30)
    return fig

def rset_classify(matrix , order):
    matrix_new = matrix[matrix[:,-1] == order[0]]
    for i in order[1:]:
        m = matrix[matrix[:,-1] == i]
        matrix_new = np.array(np.vstack((matrix_new ,m)))

    return matrix_new

def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
 

##-----------------------------get_gene_name--------------------------------------------------   
gene_type = ({'names':['gene_name'],
               'formats':['S64']})
gene_name = np.loadtxt('H:\\Workspace_New\\data\\RNA\\diff_expression\\diff_gene_matrix.txt' , dtype = gene_type , usecols = (0) , skiprows = 1)
gene_name = gene_name['gene_name']     


#------------------------------get_union_gene----------------------------------------------
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')        
               
#----------------------------------get_plot_matrix-------------------------------------------
matrix_0 = get_indexed_matrix(gene_name , False)
matrix_1 = get_indexed_matrix(gene_name , True)
matrix = K_means_cluster(matrix_0 , matrix_1[:,:-1] , 2 , True)

fig = plot_heatmap(matrix)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\gene_expression_heatmap_np.log2(FPKM+1).pdf')

outFil = open('H:\\Workspace_New\\data\\RNA\\diff_expression\\diff_gene_matrix_np.log2(FPKM+1)_classify2.txt' , 'w')
outFil.writelines('\t'.join(['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC' , 'Classify']) + '\n')
outFil0 = open('H:\\Workspace_New\\data\\RNA\\diff_expression\\diff_gene_classify2_0.bed' , 'w')
outFil1 = open('H:\\Workspace_New\\data\\RNA\\diff_expression\\diff_gene_classify2_1.bed' , 'w')
index = matrix[:,-2]
new_gene = []
n = 0
for i in index:
    clu = str(int(matrix[:,-1][n]))
    gene = union_gene[union_gene['gene_name'] == gene_name[int(i)]][0]
    outFil.writelines('\t'.join([gene['gene_name'] , gene['chr'] , gene['strand'] , str(gene['start']) , 
                                 str(gene['end']) , str(np.log2(gene['CCS'] + 1)) , str(np.log2(gene['NT5'] + 1)) , 
                                 str(np.log2(gene['NT6'] + 1)) , str(np.log2(gene['fESC'] + 1)) ,  clu]) + '\n')
    if clu == '0':
        outFil0.writelines('\t'.join([gene['chr'] , str(gene['start']) , str(gene['end']) , 'No' , str(0) , gene['strand']]) + '\n')
    else:
        outFil1.writelines('\t'.join([gene['chr'] , str(gene['start']) , str(gene['end']) , 'No' , str(0) , gene['strand']]) + '\n')
    n += 1
    
outFil.close()
outFil0.close()
outFil1.close() 
