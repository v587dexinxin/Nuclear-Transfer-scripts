# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 17:29:22 2020

@author: han-luo
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

dataFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC_ATAC\\RNA\\gene_expression'


MEF_1 = pd.read_table(os.path.join(dataFolder , 'MEF_RNA_R1.txt') , index_col = 1)
MEF_2 = pd.read_table(os.path.join(dataFolder , 'MEF_RNA_R2.txt') , index_col = 1)

MEF_1 = MEF_1[~ MEF_1.index.duplicated()]
MEF_2 = MEF_2[~ MEF_2.index.duplicated()]
union = pd.read_table('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_genes.txt' , index_col=1)

Merge = pd.concat([union[[u'Gene_ID', u'Chr', u'Start', u'End', u'CCS_1', u'CCS_2', u'CCS_3']],
           MEF_1[['FPKM']], MEF_2[['FPKM']]], axis=1)
           
Merge.dropna(subset=['Gene_ID'], inplace=True)

Merge.columns = [u'Gene_id', u'Chr', u'Start', u'End', u'CCS_1_FPKM', u'CCS_2_FPKM',
       u'CCS_3_FPKM', u'MEF_1_FPKM', u'MEF_2_FPKM']
       
Diff = pd.read_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC_ATAC\\RNA\\reads_count\\CCS_MEF.csv', index_col=1)           
           
Merged = pd.concat([Diff[[u'gene_id', u'baseMean', u'log2FoldChange', u'lfcSE', u'stat', u'pvalue', u'padj', u'MEF_is_up', u'MEF_is_down']],
                   Merge[[u'CCS_1_FPKM', u'CCS_2_FPKM', u'CCS_3_FPKM', u'MEF_1_FPKM', u'MEF_2_FPKM']]] , axis=1)


Merged = Merged.sort_values(by = ['padj'])
Merged.to_csv(os.path.join(dataFolder,'Merged_FPKM.csv'))


MEF = Merged[['MEF_1_FPKM', 'MEF_2_FPKM']].mean(axis=1)
CCS = Merged[['CCS_1_FPKM', 'CCS_2_FPKM', 'CCS_3_FPKM']].mean(axis=1)
mask =(CCS>1) | (MEF>1)

Merged = Merged[mask]
Merged.dropna(subset=['padj'], inplace=True)

index_new = []
for i in Merged.index:
    ccs = (Merged.loc[i].CCS_1_FPKM + Merged.loc[i].CCS_2_FPKM + Merged.loc[i].CCS_3_FPKM ) / 3 
    mef = (Merged.loc[i].MEF_1_FPKM + Merged.loc[i].MEF_2_FPKM ) / 2
    if ((Merged.loc[i].log2FoldChange) > 0 and (mef > ccs)) or ((Merged.loc[i].log2FoldChange) < 0 and (mef < ccs)):
        index_new.append(i)
        
 
Merged = Merged.loc[index_new]
   
Merged.to_csv('H:\\Workspace_New\\data\\IPSC_ATAC\\RNA\\reads_count\\CCS_MEF_diff_genes.csv')

Merged_select_up = Merged[(Merged.padj < 0.05) & (Merged.log2FoldChange>1.5)].index
Merged_select_down = Merged[(Merged.padj < 0.05) & (Merged.log2FoldChange<-1.5)].index
Merged_select_no = Merged[~((Merged.padj < 0.05) & (abs(Diff.log2FoldChange)>1.5))].index


plt.scatter(np.log2(CCS[Merged_select_no].values+1), np.log2(MEF[Merged_select_no].values+1), c = 'gray', alpha = 0.5)
plt.scatter(np.log2(CCS[Merged_select_up].values+1), np.log2(MEF[Merged_select_up].values+1), c = 'red', alpha = 0.5)
plt.scatter(np.log2(CCS[Merged_select_down].values+1), np.log2(MEF[Merged_select_down].values+1), c = 'blue', alpha = 0.5)
plt.xlabel('CCS')
plt.ylabel('MEF')

Diff
Diff['gene_name'] = pd.Series(Merged.index, index = Merged.Gene_id)[Diff.index]
Diff
mask_up = (Diff.padj < 0.05) & (Diff.log2FoldChange > 1.5)
mask_down = (Diff.padj < 0.05) & (Diff.log2FoldChange < -1.5)
Diff['is_up'] = mask_up
Diff['is_down'] = mask_down
Diff
Diff.to_csv('CCS_MEF_new.csv')