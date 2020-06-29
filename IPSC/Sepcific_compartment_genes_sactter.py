# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 17:44:25 2020

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
import pandas as pd

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')



    
union_gene = pd.read_csv('H:\\Workspace_New\\data\\IPSC_ATAC\\RNA\\gene_expression\\Merged_FPKM.csv' , index_col = 0)
union_gene = pd.concat([union_gene['gene_id'] , pd.DataFrame(union_gene.index , index=union_gene.index) , union_gene[['CCS_1_FPKM' , 'CCS_2_FPKM' , 'CCS_3_FPKM']].mean(1) , 
                        union_gene[['MEF_1_FPKM' , 'MEF_2_FPKM']].mean(1)] , axis = 1)
union_gene.columns = ['Gene_ID' , 'Gene_name' , 'CCs_FPKM' , 'MEF_FPKM']

union_gene_1 = pd.read_table('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt' , index_col = 0)  
 
union_gene = pd.concat([union_gene[['Gene_ID' , 'Gene_name']] , union_gene_1[[u'Chr', u'Strand', u'Start', u'End']] , union_gene[['CCs_FPKM' , 'MEF_FPKM']]] , axis = 1) 

union_gene.columns = ['Gene_ID' , 'Gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCs_FPKM' , 'MEF_FPKM']


cl1 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_A_B\\CCS_A_genes.txt' , dtype = 'S64')
cl2 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_A_B\\CCS_A_NT_A_B_genes.txt' , dtype = 'S64')
cl3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_A_B\\MEF_A_genes.txt' , dtype = 'S64')
cl4 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_A_B\\IPSC_A_B_MEF_A_genes.txt' , dtype = 'S64')

cl5 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_B_A\\CCS_B_genes.txt' , dtype = 'S64')
cl6 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_B_A\\CCS_B_NT_B_A_genes.txt' , dtype = 'S64')
cl7 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_B_A\\MEF_B_genes.txt' , dtype = 'S64')
cl8 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_B_A\\IPSC_B_A_MEF_B_genes.txt' , dtype = 'S64')



diff = pd.read_csv('H:\\Workspace_New\\data\\IPSC_ATAC\\RNA\\reads_count\\Filtered_CCS_MEF_diff_genes.csv' , index_col=0)

cl1_gene = [] ; cl2_gene = [] ; cl3_gene = [] ; cl4_gene = []
cl5_gene = [] ; cl6_gene = [] ; cl7_gene = [] ; cl8_gene = []


for i in cl1:
    if i in diff.index:
        cl1_gene.append(i)
        
        
for i in cl2:
    if i in diff.index:
        cl2_gene.append(i)
        
        
for i in cl3:
    if i in diff.index:
        cl3_gene.append(i)
        
for i in cl4:
    if i in diff.index:
        cl4_gene.append(i)

for i in cl5:
    if i in diff.index:
        cl5_gene.append(i)
        
for i in cl6:
    if i in diff.index:
        cl6_gene.append(i)

for i in cl7:
    if i in diff.index:
        cl7_gene.append(i)

for i in cl8:
    if i in diff.index:
        cl8_gene.append(i)        
        
out = pd.ExcelWriter('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_CCS_MEF_diff_genes.xlsx')
union_gene.ix[cl1_gene].to_excel(out, sheet_name='CCs_A_genes')
union_gene.ix[cl2_gene].to_excel(out, sheet_name='CCsA_NTsAB_genes')
union_gene.ix[cl3_gene].to_excel(out, sheet_name='MEF_A_genes')
union_gene.ix[cl4_gene].to_excel(out, sheet_name='MEFA_IPSCAB_genes')
union_gene.ix[cl5_gene].to_excel(out, sheet_name='CCs_B_genes')
union_gene.ix[cl6_gene].to_excel(out, sheet_name='CCsB_NTsBA_genes')
union_gene.ix[cl7_gene].to_excel(out, sheet_name='MEF_B_genes')
union_gene.ix[cl8_gene].to_excel(out, sheet_name='MEFB_IPSCBA_genes')

out.save()




pp = PdfPages('D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig4_figs\\Fig4E_Specific_A_B_genes.pdf')
size = (10, 10)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6
fig = plt.figure(figsize = size)
ax = fig.add_axes([Left , HB , width, HH])
ax.scatter(np.log2(np.array(union_gene.ix[cl1_gene].MEF_FPKM) + 1) , np.log2(np.array(union_gene.ix[cl1_gene].CCS_FPKM) + 1), alpha = 1 , c = 'blue')
ax.scatter(np.log2(np.array(union_gene.ix[cl2_gene].MEF_FPKM) + 1) , np.log2(np.array(union_gene.ix[cl2_gene].CCS_FPKM) + 1) , alpha = 1 , c = 'deeppink')
ax.scatter(np.log2(np.array(union_gene.ix[cl3_gene].MEF_FPKM) + 1) , np.log2(np.array(union_gene.ix[cl3_gene].CCS_FPKM) + 1) , alpha = 1 , c = 'lightskyblue')
ax.scatter(np.log2(np.array(union_gene.ix[cl4_gene].MEF_FPKM) + 1) , np.log2(np.array(union_gene.ix[cl4_gene].CCS_FPKM) + 1) , alpha = 1 , c = 'sandybrown')

ax.set_xlabel('MEF_log2(FPKM+1)' , size = 20)
ax.set_ylabel('CCS_log2(FPKM+1)' , size = 20)

ax.set_xlim(-0.5 , 12)
ax.set_ylim(-0.5 , 12)

ax.plot([8 , 9] , [9.5 , 9.5] , c = 'sandybrown' , lw = 2.5 )
ax.plot([8 , 9] , [10 , 10] , c = 'lightskyblue' , lw = 2.5 )
ax.plot([8 , 9] , [10.5 , 10.5] , c = 'deeppink' , lw = 2.5 )
ax.plot([8 , 9] , [11 , 11] , c = 'blue' , lw = 2.5 )

ax.text(9.3 , 9.4 , 'MEFB_IPSCBA' , size = 10 )
ax.text(9.3 , 9.9 , 'MEF_B' , size = 10 )
ax.text(9.3 , 10.4 , 'CCSB_NTBA' , size = 10 )
ax.text(9.3 , 10.9, 'CCS_B' , size = 10 )

ax.set_title('B_TO_A' , size = 20)

pp.savefig(fig)
pp.close()       