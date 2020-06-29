# -*- coding: utf-8 -*-
"""
Created on Sun Nov 24 16:57:41 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as ssd

# Our Own Color Map
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['blue' , '#FFFFFF' , 'red'])
my_cmap.set_bad('#D3D3D3')
 

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a
               
def get_union_gene(Fil):
    union_type = np.dtype({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                          'formats':['S64' , 'S64' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union_gene = list(union_gene)
    union_gene = Sort(union_gene , 1 , 3)
    union_gene = np.array(union_gene , dtype = union_type)
    return union_gene

##NT
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')  

##IPSC            
RNAFolder = 'H:\\Workspace_New\\data\\IPSC_ATAC\\RNA\\gene_expression'

cells_IPSC = {'MEF':0 , 'IPSC':1 , 'ESC':2 }
cell_IPSC = ['MEF_RNA_R1' , 'MEF_RNA_R2' , 'IPSC_RNA_R1' , 'IPSC_RNA_R2' , 'ESC_RNA_R1' , 'ESC_RNA_R2']
gene_type = np.dtype({'names':['gene_name','FPKM'],
                      'formats':['S64',np.float]})




gene_IPSC = {}
gene_name_IPSC = []
for c in cell_IPSC:
    RNAFil = c + '.txt'
    RNASource = os.path.join(RNAFolder , RNAFil)
    RNAData = np.loadtxt(RNASource , skiprows = 1 , usecols = (1 , 7) , dtype = gene_type)
    gene_IPSC[c] = RNAData
    for i in RNAData:
        if i['gene_name'] != '-':
            gene_name_IPSC.append(i['gene_name'])
        
gene_name_IPSC = set(gene_name_IPSC)

gene_all_IPSC = {}
for i in gene_name_IPSC:
    gene_all_IPSC[i] = [0 , 0 , 0]
    for c in cell_IPSC:
        ce = c.split("_")[0]
        mask = ([gene_IPSC[c]['gene_name'] == i])
        overlap = gene_IPSC[c][mask]
        if overlap.size == 0:
            fpkm = 0
        else:
            fpkm = overlap[0]['FPKM']
        gene_all_IPSC[i][cells_IPSC[ce]] += fpkm
    gene_all_IPSC[i][0] = gene_all_IPSC[i][0] / 2
    gene_all_IPSC[i][1] = gene_all_IPSC[i][1] / 2
    gene_all_IPSC[i][2] = gene_all_IPSC[i][2] / 2

    
gene_name = gene_name_IPSC.intersection(set(union_gene['gene_name']))
    
RNA_all = {'CCS':[] , 'NT5':[], 'NT6':[] , 'fESC':[] , 'MEF':[] , 'IPSC':[] , 'ESC':[]}            
for i in gene_name:
    gene_NT = union_gene[union_gene['gene_name'] == i][0]
    if (gene_NT['CCS'] != 0) or (gene_NT['NT5'] != 0) or (gene_NT['NT6'] != 0) or (gene_NT['fESC'] != 0) or (gene_all_IPSC[i][0] != 0) or (gene_all_IPSC[i][1] != 0) or (gene_all_IPSC[i][2] != 0):
        if (gene_NT['CCS'] <= 3000) and (gene_NT['NT5'] <= 3000) and (gene_NT['NT6'] <= 3000) and (gene_NT['fESC'] <= 3000) and (gene_all_IPSC[i][0] <= 3000) and (gene_all_IPSC[i][1] <= 3000) and (gene_all_IPSC[i][2] <= 3000):
            RNA_all['CCS'].append(gene_NT['CCS'])
            RNA_all['NT5'].append(gene_NT['NT5'])
            RNA_all['NT6'].append(gene_NT['NT6'])
            RNA_all['fESC'].append(gene_NT['fESC'])
            RNA_all['MEF'].append(gene_all_IPSC[i][0])
            RNA_all['IPSC'].append(gene_all_IPSC[i][1])
            RNA_all['ESC'].append(gene_all_IPSC[i][2])
            
            
            
#            RNA_all['CCS'].append(np.log2(gene_NT['CCS'] + 1))
#            RNA_all['NT5'].append(np.log2(gene_NT['NT5'] + 1))
#            RNA_all['NT6'].append(np.log2(gene_NT['NT6'] + 1))
#            RNA_all['fESC'].append(np.log2(gene_NT['fESC'] + 1))
#            RNA_all['MEF'].append(np.log2(gene_all_IPSC[i][0] + 1))
#            RNA_all['IPSC'].append(np.log2(gene_all_IPSC[i][1] + 1))
#            RNA_all['ESC'].append(np.log2(gene_all_IPSC[i][2] + 1))

    
for k , v in RNA_all.items():
    v = np.array(v)
    RNA_all[k] = v
    
  
cell = ['CCS' , 'MEF' , 'NT5' , 'NT6' , 'fESC'  , 'IPSC' , 'ESC']
cells = {'CCS':0 , 'MEF':1 , 'NT5':2 , 'NT6':3 , 'fESC':4 , 'IPSC':5 , 'ESC':6}
cor_matrix = np.zeros((7,7))
for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(RNA_all[i] , RNA_all[j])[0][1] , 3)
        


pp = PdfPages('D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S5_figs\\NT_IPSC_RNA_correlation_gene_FPKM.pdf')
left, bottom, width, height = 0.1, 0.2, 0.60, 0.6
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
im = ax.imshow(cor_matrix,vmax=0.98,vmin = 0.79,cmap=my_cmap , origin = 'lower')

x = cell
ax.set_xticks(np.arange(len(x)))
ax.set_yticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 20)
ax.set_yticklabels(x,fontsize = 20)
ax.set_title('Correlation of RNA-seq between different cells',fontsize=25)
for i in range(len(x)):
    ax.plot([-0.5 , 6.5] , [i + 0.5 , i + 0.5] , color="black")
    ax.plot([i + 0.5 , i + 0.5] , [-0.5 , 6.5] , color="black")
    for j in range(len(x)):
        text = ax.text(j, i, cor_matrix[i, j],
                       ha="center", va="center", color="black", fontsize = 20)

ax = fig.add_axes([left + width + 0.08, bottom, 0.03, height])
fig.colorbar(im, cax = ax)
pp.savefig(fig)
pp.close()


##hierarchy_cluster
data = np.vstack((RNA_all['CCS'] , RNA_all['MEF'] , RNA_all['NT5'] , RNA_all['NT6'] , RNA_all['fESC'] , RNA_all['IPSC'] , RNA_all['ESC']))
disMat = sch.distance.pdist(data,'euclidean')     
Z=sch.linkage(disMat,method='average') 
fig=sch.dendrogram(Z)