# -*- coding: utf-8 -*-
"""
Created on Wed Feb 27 20:32:54 2019

@author: xxli
"""

from __future__ import division
import numpy as np

classifies = {'CCS_noNT_nofESC':0 , 'noCCS_noNT_fESC':1 , 'noCCS_NT_fESC':2 , 'noCCS_NT_nofESC':3}
classify = ['CCS_noNT_nofESC' , 'noCCS_noNT_fESC' , 'noCCS_NT_fESC' , 'noCCS_NT_nofESC']
data_type = ({'names':['motif' , 'p-value'],
              'formats':['S64' , np.float]})

def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene

TFs = {}
TF = []
for c in classify:
    TFs[c] = []
    dataFil = 'D:\\Workspace_New\\data\\ATAC\\motif\\' + c + '\\knownResults.txt'

    data = np.loadtxt(dataFil , dtype = data_type , skiprows = 1 , usecols = (0 , 2))

    for i in data:
        motif = i['motif'].split('(')[0]
        p = i['p-value']
        if p == 0:
            l = 1000
        else:
            l = -np.log10(p)
        if p <= 1e-20:
            TFs[c].append((motif , l))
            TF.append(motif)
        else:
            continue

for k , v in TFs.items():
    v = np.array(v , dtype = data_type)
    TFs[k] = v
TF = list(set(TF))
        
union_gene = get_union_gene('D:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
p_value = np.zeros((len(TF) , 4))

for i in range(len(TF)):
    mask0 = (TFs['CCS_noNT_nofESC']['motif'] == TF[i])
    mask1 = (TFs['noCCS_noNT_fESC']['motif'] == TF[i])
    mask2 = (TFs['noCCS_NT_fESC']['motif'] == TF[i])
    mask3 = (TFs['noCCS_NT_nofESC']['motif'] == TF[i])
    overlap0 = TFs['CCS_noNT_nofESC'][mask0]
    overlap1 = TFs['noCCS_noNT_fESC'][mask1]
    overlap2 = TFs['noCCS_NT_fESC'][mask2]
    overlap3 = TFs['noCCS_NT_nofESC'][mask3]
    if  overlap0.size != 0:
        p_value[i][0] = overlap0[0]['p-value']
    else:
        p_value[i][0] = 0
        
    if  overlap1.size != 0:
        p_value[i][1] = overlap1[0]['p-value']
    else:
        p_value[i][1] = 0
    
    if  overlap2.size != 0:
        p_value[i][2] = overlap2[0]['p-value']
    else:
        p_value[i][2] = 0
        
    if  overlap3.size != 0:
        p_value[i][3] = overlap3[0]['p-value']
    else:
        p_value[i][3] = 0
    
    
    
