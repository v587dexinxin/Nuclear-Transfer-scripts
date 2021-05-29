# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 11:22:12 2021

@author: xxli
"""
'''
This code is used for the bubble chart of transcription factor binding of each type of ATAC-seq peak in fig 1F.

'''

from __future__ import division
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

classify = ['Cluster2' , 'Cluster3_Repro' , 'Cluster3_Partial' , 'Cluster3_Over']
data_type = ({'names':['motif' , 'p-value'],
              'formats':['U64' , np.float]})

def get_raw_genes_new(Fil):
    '''
    '''
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT2' , 'NT3' , 'NT4' , 'NT5' , 'NT6' , 'F35' , 'F37' , 'F40' , 'F41'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        # if (i['CCS'] > 0.1) or (i['NT5'] > 0.1) or (i['NT6'] > 0.1) or (i['F35'] > 0.1) or (i['F40'] > 0.1):
            union.append(i)
    union = np.array(union , dtype = union_type)

    return union
    
    

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

TFs = {}
for c in classify:
    TFs[c] = []
    dataFil = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\motif\\' + c + '\\knownResults.txt'

    data = np.loadtxt(dataFil , dtype = data_type , skiprows = 1 , usecols = (0 , 2))

    for i in data:
        motif = i['motif'].split('(')[0]
        p = i['p-value']
        if p == 0:
            l = 100
        else:
            if p <= 1e-100:
                l = 100
            else:              
                l = -np.log10(p)
        if p < 1: 
            TFs[c].append((motif , l))
        else:
            continue

for k , v in TFs.items():
    v = np.array(v , dtype = data_type)
    TFs[k] = v
        
union_gene = get_raw_genes_new('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_gene_expression.txt')

#Cluster2_Cluster3_Repro
gene_name = ['Mafk','Gata6','Nf1','Runx2','Runx1','Bach2','Junb' , 'Fosl2','Fosl1','Nanog','Klf9','Klf4','Pou3f1','Pou5f1','Sox3','Sox2','Ctcfl','Ctcf']
TF_name = ['MafK','Gata6','NF1','RUNX2','RUNX1','Bach2','JunB','Fra2','Fra1','Nanog','Klf9','Klf4','Oct6','Oct4','Sox3','Sox2','BORIS','CTCF']
        
p_matrix = np.zeros((len(TF_name) , 2))        
for i in range(len(TF_name)):
    mask0 = (TFs['Cluster2']['motif'] == TF_name[i])
    mask1 = (TFs['Cluster3_Repro']['motif'] == TF_name[i])
    overlap0 = TFs['Cluster2'][mask0]
    overlap1 = TFs['Cluster3_Repro'][mask1]
    if  overlap0.size != 0:
        p_matrix[i][0] = overlap0[0]['p-value']
    else:
        p_matrix[i][0] = 0
        
    if  overlap1.size != 0:
        p_matrix[i][1] = overlap1[0]['p-value']
    else:
        p_matrix[i][1] = 0

g_matrix = np.zeros((len(gene_name) , 2))        
for i in range(len(gene_name)):
    gene = union_gene[union_gene['gene_name'] == gene_name[i]]
    g_matrix[i][0] = np.log2(gene[0]['CCS'] + 1)
    g_matrix[i][1] = np.log2((gene[0]['NT5'] + gene[0]['NT6'] + gene[0]['F35'] + + gene[0]['F40'])/4 + 1)




#Bubble_Plot
sns.set(style = "whitegrid")#设置样式

x = np.array([0 , 1] * 18)
y = []
for i in range(18):
    y.extend([i,i])

z = []
for i in range(len(x)):
    z.append(p_matrix[y[i]][x[i]])#用来调整各个点的大小s
fpkm = []
for i in range(len(x)):
    fpkm.append(g_matrix[y[i]][x[i]])#用来调整各个点的大小s
    
Left = 0.3 ; HB = 0.15 ; width = 0.4 ; HH = 0.7
cm = plt.cm.get_cmap('Reds')
fig = plt.figure(figsize = (10 , 12))
ax = fig.add_axes([Left  , HB , width , HH])
bubble = ax.scatter(x, y , s = np.array(z) * 5 + 20 , c = np.array(fpkm) + 0.5, vmin = 0.5 , vmax = 8.5 , cmap = cm, linewidth = 0.5 , edgecolors = 'black')
#ax.scatter([x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5, alpha = 0.8)
ax.grid()
    
ax.set_xticks([-1, 0 , 1 , 2 ])
ax.set_xticklabels(['' , 'CCS' , 'NT_fESC' , ''] , fontsize = 15)
ax.set_yticks(range(18))
ax.set_yticklabels(gene_name , fontsize = 20)
ax.set_xlabel('classifies', fontsize = 25 , labelpad = 40)#X轴标签
ax.set_ylabel('TFs', fontsize = 25)#Y轴标签

ax2 = fig.add_axes([Left + 0.34 , HB - 0.05 , 0.06 , 0.015])
fig.colorbar(bubble,cax = ax2, orientation='horizontal' , ticks = [0.5, 8.5])

    
    
ax = fig.add_axes([Left + width + 0.05  , HB , 0.2 , 0.15])
ax.scatter([1,1,1] , [1,1.45,2.0] ,s = [20 , 270 , 520] , marker = 'o' , c = 'white' ,linewidth = 0.5,edgecolors = 'black')
ax.set_ylim((0.7,2.3))
ax.set_yticks([1 , 1.5 ,2])
ax.set_yticklabels(['0' , '50' , '100'])
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\motif\\bubble_plot\\Cluster2_Cluster3_Repro_motif_bubble.pdf')


#######################################################
######Cluster3_Resis,Cluster3_Over

classify = ['Cluster3_Partial','Cluster3_Over']
data_type = ({'names':['motif' , 'p-value'],
              'formats':['U64' , np.float]})
TFs = {}
for c in classify:
    TFs[c] = []
    dataFil = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\motif\\' + c + '\\knownResults.txt'

    data = np.loadtxt(dataFil , dtype = data_type , skiprows = 1 , usecols = (0 , 2))

    for i in data:
        motif = i['motif'].split('(')[0]
        p = i['p-value']
        if p == 0:
            l = 500
        else:
            l = -np.log10(p)
        if p < 1: 
            TFs[c].append((motif , l))
        else:
            continue

for k , v in TFs.items():
    v = np.array(v , dtype = data_type)
    TFs[k] = v
        







gene_name = ['Hinfp','Tcf3','Zbtb12','Snai2','Snai1','Trp53','Yy1','Zeb2','Zeb1','Zic3','Zic2','Rilp']
TF_name = ['HINFP','E2A','ZBTB12','Slug','Snail1','p53','YY1','ZEB2','ZEB1','Zic3','Zic2','REST-NRSF']
 
p_matrix = np.zeros((len(TF_name) , 2))        
for i in range(len(TF_name)):
    mask0 = (TFs['Cluster3_Partial']['motif'] == TF_name[i])
    mask1 = (TFs['Cluster3_Over']['motif'] == TF_name[i])
    overlap0 = TFs['Cluster3_Partial'][mask0]
    overlap1 = TFs['Cluster3_Over'][mask1]
    if  overlap0.size != 0:
        p_matrix[i][0] = overlap0[0]['p-value']
    else:
        p_matrix[i][0] = 0
        
    if  overlap1.size != 0:
        p_matrix[i][1] = overlap1[0]['p-value']
    else:
        p_matrix[i][1] = 0

g_matrix = np.zeros((len(gene_name) , 2))        
for i in range(len(gene_name)):
    gene = union_gene[union_gene['gene_name'] == gene_name[i]]    
    g_matrix[i][0] = np.log2((gene[0]['F35'] + gene[0]['F40']) /2 + 1)
    g_matrix[i][1] = np.log2((gene[0]['NT5'] + gene[0]['NT6']) /2 + 1)




#Bubble_Plot
sns.set(style = "whitegrid")#设置样式

x = [0,1] * 12
y = []
for i in range(12):
    y.extend([i,i])

z = []
for i in range(len(x)):
    if p_matrix[y[i]][x[i]] > 50:
        z.append(50)
    else:
        z.append(p_matrix[y[i]][x[i]])#用来调整各个点的大小s
fpkm = []
for i in range(len(x)):
    fpkm.append(g_matrix[y[i]][x[i]])#用来调整各个点的大小s
    
Left = 0.3 ; HB = 0.15 ; width = 0.4 ; HH = 0.4
cm = plt.cm.get_cmap('Reds')
fig = plt.figure(figsize = (10 , 12))
ax = fig.add_axes([Left  , HB , width , HH])
bubble = ax.scatter(x, y , s = np.array(z) * 10 + 20, c = np.array(fpkm) + 0.5, vmin = 0.5 , vmax = 8.5 , cmap = cm, linewidth = 0.5 , edgecolors = 'black')
#ax.scatter([x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5, alpha = 0.8)
ax.grid()
    
ax.set_xticks([-1, 0 , 1 , 2 ])
ax.set_xticklabels(['' , 'Cluster3_Partial' , 'Cluster3_Over' , ''] , fontsize = 15)
ax.set_yticks(range(12))
ax.set_yticklabels(gene_name , fontsize = 20)
ax.set_xlabel('classifies', fontsize = 25 , labelpad = 40)#X轴标签
ax.set_ylabel('TFs', fontsize = 25)#Y轴标签

ax2 = fig.add_axes([Left + 0.34 , HB - 0.05 , 0.06 , 0.015])
fig.colorbar(bubble,cax = ax2, orientation='horizontal' , ticks = [0.5 , 8.5])

    
    
ax = fig.add_axes([Left + width + 0.05  , HB , 0.2 , 0.15])
ax.scatter([1,1,1] , [1.0, 1.45,2.0] ,s = [20 , 270 , 520] , marker = 'o' , c = 'white' ,linewidth = 0.5,edgecolors = 'black')
ax.set_ylim((0.7,2.3))
ax.set_yticks([1 , 1.5 ,2])
ax.set_yticklabels(['0' , '25' , '50'])
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\motif\\bubble_plot\\Cluster3_Resis_Cluster3_Over_motif_bubble.pdf')

