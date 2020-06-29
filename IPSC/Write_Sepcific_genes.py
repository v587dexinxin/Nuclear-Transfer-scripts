# -*- coding: utf-8 -*-
"""
Created on Thu May 07 13:55:47 2020

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

union_gene = pd.concat([union_gene['gene_id'] , pd.DataFrame(union_gene.index , index=union_gene.index)] , axis = 1)
union_gene.columns = ['Gene_ID' , 'Gene_name']

union_gene_1 = pd.read_table('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt' , index_col = 0)  
 
union_gene = pd.concat([union_gene[['Gene_ID' , 'Gene_name']] , union_gene_1[[u'Chr', u'Strand', u'Start', u'End']]] , axis = 1) 

union_gene.columns = ['Gene_ID' , 'Gene_name' , 'chr' , 'strand' , 'start' , 'end']


cl1 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_A_B\\NT_A_B_genes.txt' , dtype = 'S64')
cl2 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_A_B\\IPSC_A_B_genes.txt' , dtype = 'S64')
cl3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_B_A\\NT_B_A_genes.txt' , dtype = 'S64')
cl4 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_B_A\\IPSC_B_A_genes.txt' , dtype = 'S64')

cl5 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_resist_A_B\\Resist_NT_A_B_genes.txt' , dtype = 'S64')
cl6 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_resist_A_B\\Resist_IPSC_A_B_genes.txt' , dtype = 'S64')
cl7 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_resist_B_A\\Resist_NT_B_A_genes.txt' , dtype = 'S64')
cl8 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\Specific_resist_B_A\\Resist_IPSC_B_A_genes.txt' , dtype = 'S64')


  
        
out = pd.ExcelWriter('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\chromatin_switch_&_resist_specific_genes.xlsx')
union_gene.ix[cl1].to_excel(out, sheet_name='NTs_specificAB_genes')
union_gene.ix[cl2].to_excel(out, sheet_name='IPSC_specificAB_genes')
union_gene.ix[cl3].to_excel(out, sheet_name='NTs_specificBA_genes')
union_gene.ix[cl4].to_excel(out, sheet_name='IPSC_specificBA_genes')
union_gene.ix[cl5].to_excel(out, sheet_name='Resis_NTs_specificAB_genes')
union_gene.ix[cl6].to_excel(out, sheet_name='Resis_IPSC_specificAB_genes')
union_gene.ix[cl7].to_excel(out, sheet_name='Resis_NTs_specificBA_genes')
union_gene.ix[cl8].to_excel(out, sheet_name='Resis_IPSC_specificBA_genes')

out.save()

