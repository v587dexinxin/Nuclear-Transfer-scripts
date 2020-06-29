# -*- coding: utf-8 -*-
"""
Created on Thu Dec 06 19:28:25 2018

@author: xxli
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
                                            ['skyblue','k' , 'yellow'])
my_cmap.set_bad('#2672a1')
                

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

def K_means_cluster(matrix_0,matrix,n):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1]))
    return matrix_1

    
def plot_heatmap(matrix):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    vmax = np.percentile(matrix,95)
    vmin = np.percentile(matrix,5)
    im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    plt.colorbar(im)
    x = ['CCS','NT5','NT6','fESC']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 30)
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



    
CCS_NT = get_common2_gene_name('D:/Workspace_New/data/RNA/diff_expression/basemean_500_fc_1.5/Filtered_CCS_NT5.txt' , 
                               'D:/Workspace_New/data/RNA/diff_expression/basemean_500_fc_1.5/Filtered_CCS_NT6.txt')
fESC_NT = get_common2_gene_name('D:/Workspace_New/data/RNA/diff_expression/basemean_500_fc_1.5/Filtered_NT5_fESC.txt' , 
                               'D:/Workspace_New/data/RNA/diff_expression/basemean_500_fc_1.5/Filtered_NT6_fESC.txt')
CCS_fESC = get_union_gene_name(['CCS_fESC'])


union_gene = get_union_gene('D:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

gene_name = list(CCS_NT.union(fESC_NT).union(set(CCS_fESC)))

#c1 = []
#c2 = []
#c3 = []
# 
#for i in CCS_fESC:
#    if i not in fESC_NT:
#        if i in CCS_NT:
#            c2.append(i)
#        else:
#            c1.append(i)
#            
#for i in CCS_NT:
#    if i not in fESC_NT:
#        if i not in CCS_fESC:
#            c3.append(i)
#        else:
#            continue

#al = {'CCS_fESC':c1 , 'CCS_fESC_NT':c2 , 'CCS_NT':c3}
#
#cells = ['CCS_fESC_+' , 'CCS_fESC_NT_+' , 'CCS_NT_+' , 'CCS_fESC_-' , 'CCS_fESC_NT_-' , 'CCS_NT_-']
#gene_names = {}
#for c in cells:
#    gene_names[c] = []
#for i in al:
#    for j in al[i]: 
#        gene = union_gene[union_gene['gene_name'] == j][0]
#        if gene['CCS'] > gene['fESC']:
#            gene_names[i + '_+'].append(j)
#        else:
#            gene_names[i + '_-'].append(j)
#    
#gene_name = []
#for c in cells:
#    for j in gene_names[c]:
#        gene_name.append(j)
                 
matrix = get_indexed_matrix(gene_name , zscore=True)
matrix_new = K_means_cluster(matrix , matrix[:,:-1] , 10)
rmatrix = rset_classify(matrix_new ,order=[6,0,8,4,3,7,1,9,2,5])
plot_heatmap(rmatrix[:,:-2])

gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')

o = open('D:\\Workspace_New\\data\\RNA\\diff_expression\\diff_gene_matrix.txt','w')
o.writelines('\t'.join(['Gene_name' , 'Chr' , 'Start' , 'End' , 'Strand' , 'CCS' , 'NT5' 'NT6' , 'fESC' , 'Classify']) + '\n')
for i in rmatrix:
    g_name = gene_name[int(i[4])]
    classify = int(i[-1])
    gene = gtf[gtf['gene_name'] == g_name]
    if (classify == 6) or (classify == 0) or (classify == 8):
        c = 0
    elif(classify == 4):
        c = 1
    elif(classify == 3):
        c = 2
    elif(classify == 7) or (classify == 1):
        c = 3
    elif(classify == 9) or (classify == 2):
        c = 4
    elif(classify == 5):
        c = 5
    o.writelines('\t'.join([gene[0]['gene_name'] , gene[0]['chr'] , str(gene[0]['start']) , str(gene[0]['end']) , gene[0]['strand'] , 
                            str(i[0]) , str(i[1]) , str(i[2]) , str(i[3]) , str(c)]) + '\n')
o.close()
        