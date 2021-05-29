# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 16:14:40 2020

@author: xxli
"""


import numpy as np
import pandas as pd


def get_genes(Fil):
    '''
    '''
    data = pd.read_table(Fil , index_col = 0 , header=0)
    
    return data

def get_union_gene_new_ave(union_gene_new):
    '''
    '''
    for c in['F35' , 'F37' , 'F41' , 'NT2' , 'NT3' , 'NT4']:
        tmp = union_gene_new[[c + "_R1" , c + "_R2" , c + "_R3"]]
        union_gene_new[c + "_FPKM"] = tmp.mean(axis=1)
        
    return union_gene_new
        
    
    
    
union_gene_raw = get_genes('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')
union_gene_new = get_genes('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_additional_experiments\\FPKM\\RNA_new_gene_expression_all.txt')
union_gene_ave = get_union_gene_new_ave(union_gene_new)

union_gene = pd.concat([union_gene_raw , union_gene_ave[['NT2_FPKM','NT3_FPKM','NT4_FPKM','F35_FPKM','F37_FPKM','F41_FPKM']]] , axis = 1)

union_gene = union_gene.reindex(columns=['Chr' , 'Strand' , 'Start' , 'End' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' ,                                         
                                         'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'fESC_FPKM' , 'F41_FPKM' ])

union_gene.columns=['Chr' , 'Strand' , 'Start' , 'End' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' ,                                         
                                         'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM' ]


union_gene.to_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_gene_expression.csv' , sep=',')

















