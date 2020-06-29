# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 14:36:30 2019

@author: han-luo
"""

    
from __future__ import division
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

classify = ['union_cluster2' , 'union_cluster3']
data_type = ({'names':['motif' , 'p-value'],
              'formats':['S64' , np.float]})

def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
    

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

TFs = {}
for c in classify:
    TFs[c] = []
    dataFil = 'H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\motif\\' + c + '\\knownResults.txt'

    data = np.loadtxt(dataFil , dtype = data_type , skiprows = 1 , usecols = (0 , 2))

    for i in data:
        motif = i['motif'].split('(')[0]
        p = i['p-value']
        if p == 0:
            l = 1000
        else:
            l = -np.log10(p)
        if p < 1: 
            TFs[c].append((motif , l))
        else:
            continue

for k , v in TFs.items():
    v = np.array(v , dtype = data_type)
    TFs[k] = v
        
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

#Cluster2_cluster3
gene_name = ['Batf','Atf2','Gata1','Bach2','Fosl2','Fosl1','Runx2','Runx1','Junb','Jun','Tcf3','Zeb1','Nanog','Zic3','Pou3f1','Pou5f1','Sox3','Sox2','Klf4','Klf6']
TF_name = ['BATF','Atf2','Gata1','Bach2','Fosl2','Fra1','RUNX2','RUNX1','JunB','Jun','E2A','ZEB1','Nanog','Zic3','Oct6','Oct4','Sox3','Sox2','Klf4','KLF6']
        
p_matrix = np.zeros((len(TF_name) , 2))        
for i in range(len(TF_name)):
    mask0 = (TFs['union_cluster2']['motif'] == TF_name[i])
    mask1 = (TFs['union_cluster3']['motif'] == TF_name[i])
    overlap0 = TFs['union_cluster2'][mask0]
    overlap1 = TFs['union_cluster3'][mask1]
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
    g_matrix[i][1] = np.log2((gene[0]['NT5'] + gene[0]['NT6'] + gene[0]['fESC'])/3 + 1)




#Bubble_Plot
sns.set(style = "whitegrid")#设置样式

x = range(2) * 20
y = []
for i in range(20):
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
bubble = ax.scatter(x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5 , edgecolors = 'black')
#ax.scatter([x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5, alpha = 0.8)
ax.grid()
    
ax.set_xticks([-1, 0 , 1 , 2 ])
ax.set_xticklabels(['' , 'CCS' , 'NT_fESC' , ''] , fontsize = 15)
ax.set_yticks(range(20))
ax.set_yticklabels(gene_name , fontsize = 20)
ax.set_xlabel('classifies', fontsize = 25 , labelpad = 40)#X轴标签
ax.set_ylabel('TFs', fontsize = 25)#Y轴标签

ax2 = fig.add_axes([Left + 0.34 , HB - 0.05 , 0.06 , 0.015])
fig.colorbar(bubble,cax = ax2, orientation='horizontal' , ticks = [np.round(min(fpkm),1) + 0.6, int(max(fpkm)/2) , int(max(fpkm))])

    
    
ax = fig.add_axes([Left + width + 0.05  , HB , 0.2 , 0.15])
ax.scatter([1,1,1] , [1,1.45,2.0] ,s = [5 , 255 , 505] , marker = 'o' , c = 'white' ,linewidth = 0.5,edgecolors = 'black')
ax.set_ylim((0.7,2.3))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\Allpeaks_classify_motif_bubble.pdf')


#######################################################
######cluster3b_cluster3c

classify = ['union_cluster3_c2' , 'union_cluster3_c3']
data_type = ({'names':['motif' , 'p-value'],
              'formats':['S64' , np.float]})
TFs = {}
for c in classify:
    TFs[c] = []
    dataFil = 'H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\union_new\\cluster1_4_23\\cluster3\\motif\\' + c + '\\knownResults.txt'

    data = np.loadtxt(dataFil , dtype = data_type , skiprows = 1 , usecols = (0 , 2))

    for i in data:
        motif = i['motif'].split('(')[0]
        p = i['p-value']
        if p == 0:
            l = 1000
        else:
            l = -np.log10(p)
        if p < 1: 
            TFs[c].append((motif , l))
        else:
            continue

for k , v in TFs.items():
    v = np.array(v , dtype = data_type)
    TFs[k] = v
        
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

genenames = []
for c in TFs:
    for i in TFs[c]:
        if i['motif'] not in union_gene['gene_name']:
            print i             
        else:
            genenames.append(i['motif'])

'''
genenames.append(...)
'''
TFnames = []
for i in genenames:
    if (i in TFs['union_cluster3_c2']['motif']) or (i in TFs['union_cluster3_c3']['motif']):
        TFnames.append(i)
    else:
        print i

'''
TFnames.append(...)

TFnames : 'Foxo3,Foxf1,Eomes,Max,Ascl1,Ptf1a,Tbx5,Tgif2,Usf2,Smad3,Meis1,Cdx2,Pax7,Barx1,Zic3,Sox3,Sox10,Sox2,Sox15,Sox6,Sox17,Sox4,Fli1,Sox9,Etv2,Otx2,Hand2,Elk1,Smad2,Tcf3,Nanog,Smad4,Bcl6,Elk4,Tcf4,Tgif1,Foxo1,Tbx5,Pknox1,Tcf12,Foxa3,Pbx3,Gfi1b,Rfx2,Gata6,Rfx1,Hoxb4,ZBTB12,Slug,NRF1,FOXP1,MafK,MyoG,CTCF,p63,MafB,BORIS,NPAS2,Srebp1a,Srebp2,THRb,E2A,E2F3,AP-2alpha,ERG,ETS1,GABPA,ZEB1,ETV1,EHF,ELF3,SPDEF,YY1,p53,p73'
genenames :'Foxo3,Foxf1,Eomes,Max,Ascl1,Ptf1a,Tbx5,Tgif2,Usf2,Smad3,Meis1,Cdx2,Pax7,Barx1,Zic3,Sox3,Sox10,Sox2,Sox15,Sox6,Sox17,Sox4,Fli1,Sox9,Etv2,Otx2,Hand2,Elk1,Smad2,Tcf3,Nanog,Smad4,Bcl6,Elk4,Tcf4,Tgif1,Foxo1,Tbx5,Pknox1,Tcf12,Foxa3,Pbx3,Gfi1b,Rfx2,Gata6,Rfx1,Hoxb4,Zbtb12,Snai2,Nrf1,Foxp1,Mafk,Myog,Ctcf,Trp63,Mafb,Ctcfl,Npas2,Srebf1,Srebf2,Thrb,Tcf3,E2f3,Tfap2a,Erg,Ets1,Gabpa,Zeb1,Etv1,Ehf,Elf3,Spdef,Yy1,Trp53,Trp73'
'''

gene_name = []
TF_name = []
for i in range(len(genenames)):
    gene = union_gene[union_gene['gene_name'] == genenames[i]][0]
    TF = TFnames[i]
    c2 = TFs['union_cluster3_c2'][TFs['union_cluster3_c2']['motif'] == TF]
    c3 = TFs['union_cluster3_c3'][TFs['union_cluster3_c3']['motif'] == TF]
    if (c2.size != 0) and (c3.size != 0):
        if c2[0]['p-value'] > c3[0]['p-value']:
            if (gene['CCS'] < gene['NT5']) and (gene['CCS'] < gene['NT6']) and (gene['CCS'] < gene['fESC']):
                if (gene['fESC'] < gene['NT5']) or (gene['fESC'] < gene['NT6']):
                    gene_name.append(genenames[i])
                    TF_name.append(TFnames[i])
        elif c2[0]['p-value'] < c3[0]['p-value']:
            if (gene['CCS'] < gene['NT5']) and (gene['CCS'] < gene['NT6']) and (gene['CCS'] < gene['fESC']):            
                if (gene['fESC'] > gene['NT5']) or (gene['fESC'] > gene['NT6']):
                    gene_name.append(genenames[i])
                    TF_name.append(TFnames[i])
        else:
            pass
    elif (c2.size != 0) and (c3.size == 0):
        if (gene['CCS'] < gene['NT5']) and (gene['CCS'] < gene['NT6']) and (gene['CCS'] < gene['fESC']):            
                if (gene['fESC'] < gene['NT5']) or (gene['fESC'] < gene['NT6']):
                    gene_name.append(genenames[i])
                    TF_name.append(TFnames[i])
    elif (c2.size == 0) and (c3.size != 0):
        if (gene['CCS'] < gene['NT5']) and (gene['CCS'] < gene['NT6']) and (gene['CCS'] < gene['fESC']):
                if (gene['fESC'] > gene['NT5']) or (gene['fESC'] > gene['NT6']):
                    gene_name.append(genenames[i])
                    TF_name.append(TFnames[i])
    else:
        pass
            
            
TF_name_0 = ['YY1' , 'p53','Slug','NRF1','CTCF','Eomes','Zic3','Sox3','Sox2','Sox15','Tcf3','Nanog','ERG','GABPA']            
gene_name_0 = ['Yy1' , 'Trp53','Snai2','Nrf1','Ctcf','Eomes','Zic3','Sox3','Sox2','Sox15','Erg','Gabpa']        
gene_name = []
TF_name = []
for i in range(len(gene_name_0) - 1 , -1 , -1):
    gene_name.append(gene_name_0[i])
    TF_name.append(TF_name_0[i])

p_matrix = np.zeros((len(TF_name) , 2))        
for i in range(len(TF_name)):
    mask0 = (TFs['union_cluster3_c2']['motif'] == TF_name[i])
    mask1 = (TFs['union_cluster3_c3']['motif'] == TF_name[i])
    overlap0 = TFs['union_cluster3_c2'][mask0]
    overlap1 = TFs['union_cluster3_c3'][mask1]
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
    g_matrix[i][0] = np.log2((gene[0]['NT5'] + gene[0]['NT6']) /2 + 1)
    g_matrix[i][1] = np.log2(gene[0]['fESC'] + 1)




#Bubble_Plot
sns.set(style = "whitegrid")#设置样式

x = range(2) * 12
y = []
for i in range(12):
    y.extend([i,i])

z = []
for i in range(len(x)):
    z.append(p_matrix[y[i]][x[i]])#用来调整各个点的大小s
fpkm = []
for i in range(len(x)):
    fpkm.append(g_matrix[y[i]][x[i]])#用来调整各个点的大小s
    
Left = 0.3 ; HB = 0.15 ; width = 0.4 ; HH = 0.4
cm = plt.cm.get_cmap('Reds')
fig = plt.figure(figsize = (10 , 12))
ax = fig.add_axes([Left  , HB , width , HH])
bubble = ax.scatter(x, y , s = np.array(z) * 10 + 5, c = np.array(fpkm) + 0.5, vmin = 0.6 , vmax = 9.0 , cmap = cm, linewidth = 0.5 , edgecolors = 'black')
#ax.scatter([x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5, alpha = 0.8)
ax.grid()
    
ax.set_xticks([-1, 0 , 1 , 2 ])
ax.set_xticklabels(['' , 'cluster3_c2' , 'cluster3_c3' , ''] , fontsize = 15)
ax.set_yticks(range(12))
ax.set_yticklabels(gene_name , fontsize = 20)
ax.set_xlabel('classifies', fontsize = 25 , labelpad = 40)#X轴标签
ax.set_ylabel('TFs', fontsize = 25)#Y轴标签

ax2 = fig.add_axes([Left + 0.34 , HB - 0.05 , 0.06 , 0.015])
fig.colorbar(bubble,cax = ax2, orientation='horizontal' , ticks = [0.6 , 9.0])

    
    
ax = fig.add_axes([Left + width + 0.05  , HB , 0.2 , 0.15])
ax.scatter([1,1,1] , [1.0, 1.45,2.0] ,s = [5 , 255 , 505] , marker = 'o' , c = 'white' ,linewidth = 0.5,edgecolors = 'black')
ax.set_ylim((0.7,2.3))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\Allpeaks_cluster3_motif_bubble.pdf')

