# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:33:08 2019

@author: han-luo
"""

from __future__ import division
import math
import numpy as np
import pandas as pd
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


def get_union_gene(Fil):
    union_type = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                          'formats':['S64' , 'S64' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union_gene = list(union_gene)
    union_gene = Sort(union_gene , 1 , 3)
    union_gene = np.array(union_gene , dtype = union_type)
    return union_gene
    
def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_name = i.strip().split('gene_name')[1].split('\"')[1]
            gene_id = i.strip().split('gene_id')[1].split('\"')[1]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf.close()
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf
    
    
def get_intergenic(union_gene):
    inter_type = np.dtype({'names':['chr' , 'start' , 'end' , 'gene1' , 'gene2' ],
                           'formats':['S64' , np.int , np.int , 'S64' , 'S64']})

    intergene = []
    chrom = set(union_gene['chr'])
    chroms = []
    for g in chrom:
        if 'chr' in g:
            chroms.append(g)
    chroms.remove('chrM')
    for g in chroms:
        genes = union_gene[union_gene['chr'] == g]
        intergene.append((g , 0 , genes[0]['start'] , genes[0]['gene_name'] , genes[0]['gene_name']))
        for i in range(len(genes) - 1):
            inter = (g , genes[i]['end'] , genes[i + 1]['start'] , genes[i]['gene_name'] , genes[i + 1]['gene_name'])
            intergene.append(inter)
        intergene.append((g , genes[i + 1]['end'] , mm10[mm10['chr'] == g][0]['length'] , genes[i + 1]['gene_name'] , genes[i + 1]['gene_name']))
    intergene = np.array(intergene , dtype = inter_type)
    return intergene
            
    
def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a

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
    numbers = {'1':len(CCS_gene) , '2':len(ESC_gene)}

    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 30)
    ax.set_ylabel('\n'.join(['cluster' + str(int(i)) + ':' + str(numbers[str(int(i))]) for i in [1,2]]) , fontsize = 20)
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('Zscore ,gene numbers:' + str(len(matrix)),fontsize = 30)
    return fig

    
def get_peaks_nearest_gene(peaks , union_gene , intergenic):
    Genes = []
    for i in peaks:
        chro = i['chr']
        start = i['start']
        end = i['end']
        genes = union_gene[union_gene['chr'] == chro]
        inters = intergenic[intergenic['chr'] == chro]
        mask = (genes['start'] <= end) & (genes['end'] >= start)
        overlap = genes[mask]
        if overlap.size == 1:
            gene = overlap[0]
        elif overlap.size > 1:
            dist = {}
            for j in overlap:
                d = min(abs(j['start'] - start) , abs(j['end'] - start) , abs(j['start'] - end) , abs(j['end'] - end))
                dist[d] = j
            gene = dist[min(dist.keys())]
        else:
            mask1 = (inters['start'] <= start) & (inters['end'] >= end)
            overlap1 = inters[mask1]
            if overlap1.size == 0:
                print i 
            d1 = abs(start - overlap1[0]['start'])
            d2 = abs(overlap1[0]['end'] - end)
            if d1 <= d2:
                gene = union_gene[union_gene['gene_name'] == overlap1[0]['gene1']][0]
            else:
                gene = union_gene[union_gene['gene_name'] == overlap1[0]['gene2']][0]
        Genes.append(gene)
    Genes = np.array(Genes , dtype = union_gene.dtype)
        
    return Genes            
        
peak_type = np.dtype({'names':['chr' , 'start' , 'end'],
                      'formats':['S64' , np.int , np.int]})
lenth_type = np.dtype({'names':['chr' , 'length'],
                      'formats':['S64' , np.int ]})
union_type = np.dtype({'names':['gene_name' , 'chr' , 'starnd' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                       'formats':['S64' , 'S64' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})


def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

#get_genom                      
mm10 = np.loadtxt('G:\\data\\genome\\mm10_chr.txt' , dtype = lenth_type)
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')  
intergenic = get_intergenic(union_gene)
gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')

#get_peaks
cluster2 = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\union_cluster2.bed' , dtype = peak_type)
cluster3 = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\union_cluster3.bed' , dtype = peak_type)

#get_peak_nearest_gene
cluster2_genes = get_peaks_nearest_gene(cluster2 , union_gene , intergenic)
cluster3_genes = get_peaks_nearest_gene(cluster3 , union_gene , intergenic)


#get_diff_genes
CCS_NT5 = pd.read_csv('H:\\Workspace_New\\data\\RNA\\diff_expression\\New_fc_1.5\\Filtered_CCS_NT5.csv')
CCS_NT6 = pd.read_csv('H:\\Workspace_New\\data\\RNA\\diff_expression\\New_fc_1.5\\Filtered_CCS_NT6.csv')
CCS_fESC = pd.read_csv('H:\\Workspace_New\\data\\RNA\\diff_expression\\New_fc_1.5\\Filtered_CCS_fESC.csv')
diff_gene_names = []
for i in [CCS_NT5 , CCS_NT6 , CCS_fESC]:
    for j in i.gene_name:
        diff_gene_names.append(j)
        
diff_gene_names = set(diff_gene_names)

CCS_gene = []
for i in cluster2_genes:
    if (i['CCS'] / i['NT5'] > 1.5) and (i['CCS'] / i['NT6'] > 1.5) and (i['CCS'] / i['fESC'] > 1.5):
        if i['gene_name'] in diff_gene_names:
            CCS_gene.append(i)
CCS_gene = np.array(CCS_gene , dtype = union_type)
        
ESC_gene = []
for i in cluster3_genes:
    if (i['NT5'] / i['CCS'] > 1.5) and (i['NT6'] / i['CCS'] > 1.5) and (i['fESC'] / i['CCS'] > 1.5):
        if i['gene_name'] in diff_gene_names:
            ESC_gene.append(i)
ESC_gene = np.array(ESC_gene , dtype = union_type)


#remove replicates
CCS_gene_names = list(set(CCS_gene['gene_name']))
ESC_gene_names = list(set(ESC_gene['gene_name']))

CCS_gene_id = []
for i in CCS_gene_names:
    g_id = gtf[gtf['gene_name'] == i][0]['gene_id']
    CCS_gene_id.append(g_id)
 
ESC_gene_id = []
for i in ESC_gene_names:
    g_id = gtf[gtf['gene_name'] == i][0]['gene_id']
    ESC_gene_id.append(g_id) 

w = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\diff_genes\\union_cluster2_geneID.txt' , 'w')
for i in CCS_gene_id:    
    w.writelines(i + '\n')
w.close()


w = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\diff_genes\\union_cluster3_gene.txt' , 'w')
for i in ESC_gene_names:
    w.writelines(i + '\n')
w.close()   


#re_get_geneDatas
CCS_gene = []
for i in CCS_gene_names:
    gene = union_gene[union_gene['gene_name'] == i][0]
    CCS_gene.append(gene)
CCS_gene = np.array(CCS_gene , dtype = union_type)   

ESC_gene = []
for i in ESC_gene_names:
    gene = union_gene[union_gene['gene_name'] == i][0]
    ESC_gene.append(gene)
ESC_gene = np.array(ESC_gene , dtype = union_type)   


  
##Get_FPKM_Matrix
matrix = np.zeros((len(CCS_gene) + len(ESC_gene) , 4))

for i in range(len(CCS_gene_names)):
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
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\Selected_diff_gene_heatmap_all_peaks_Zscore_bwr.pdf')
    
    

CCS_gene = pd.DataFrame(CCS_gene)
CCS_gene = CCS_gene[['gene_name' , 'CCS' , 'NT5' , 'NT6' , 'fESC']]

ESC_gene = pd.DataFrame(ESC_gene)
ESC_gene = ESC_gene[['gene_name' , 'CCS' , 'NT5' , 'NT6' , 'fESC']]

writer = pd.ExcelWriter('D:\\ntESC_3Dreprogramming\\supp_Tables\\supp_Table3.xlsx')
CCS_gene.to_excel(writer, sheet_name='Cluster2_genes')
ESC_gene.to_excel(writer, sheet_name='Cluster3_genes')
NT_gene.to_excel(writer, sheet_name='Cluster3b_genes')
ESC_gene_1.to_excel(writer, sheet_name='Cluster3c_genes')
writer.save()



    