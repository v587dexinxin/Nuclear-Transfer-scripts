# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 21:42:42 2020

@author: xxli
"""


import pandas as pd


def Get_FPKM(df , cell_list):
    gene_names = list(df.index)
    FPKM = union_gene.loc[gene_names , cell_list]
    data = pd.concat([df, FPKM] , axis = 1)
    return data

    

df1 = pd.read_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\diff_expression\\All_RNA_ESC_vs_NTs_diff_expression_q_0.05_fc1.5.csv' , index_col = 1 , header = 0)
df2 = pd.read_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\diff_expression\\ESC_vs_NTs_diff_expression_q_0.05_fc1.5.csv' , index_col = 1 , header = 0)
union_gene = pd.read_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\gene_expression\\all_gene_expression.csv' , index_col = 0 , header = 0)

union_gene = union_gene[~union_gene.index.duplicated(keep='first')]

n = 0
for i in df1.index:
    if i in df2.index:
        n += 1
        print(i)
        
        
        

data1 = Get_FPKM(df1 , ['CCS_FPKM', 'F35_FPKM', 'F37_FPKM', 'F40_FPKM', 'F41_FPKM', 'NT2_FPKM', 'NT3_FPKM', 'NT4_FPKM', 'NT5_FPKM', 'NT6_FPKM'])
data1.to_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\diff_expression\\All_RNA_ESC_vs_NTs_diff_expression_q_0.05_fc1.5_1.csv'  , sep=',')

data2 = Get_FPKM(df2 , ['CCS_FPKM' , 'F35_FPKM' , 'F40_FPKM' , 'NT5_FPKM' , 'NT6_FPKM'])
data2.to_csv('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\RNA\\RNA_New\\diff_expression\\ESC_vs_NTs_diff_expression_q_0.05_fc1.5_1.csv'  , sep=',')






