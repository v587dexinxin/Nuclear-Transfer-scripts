# -*- coding: utf-8 -*-
"""
Created on Thu Dec 06 17:00:20 2018

@author: xxli
"""

import csv
import numpy as np
from itertools import islice

cell = ['CCS_fESC' , 'CCS_NT5' , 'CCS_NT6' , 'fESC_NT5' , 'fESC_NT6' , 'NT5_NT6']

for c in cell:
    c1 = c.split("_")[0] ; c2 = c.split("_")[1]
    gene_type = ({'names':['gene_id' , 'baseMean' , 'log2FoldChange' , 'lfcSE' , 'stat' , 'pvalue' , 'padj' , 'gene_name' , 'chr' , 'strand' , 
                      'start' , 'end' , c1  , c2 ],
                 'formats':['S64' , np.float , np.float , np.float , np.float , np.float , np.float , 'S64' , 'S8' , 'S8' , np.int , np.int , 
                            np.float , np.float]})
    f = open('F:\\xxli\\data_analysis\\BDF1\\RNA_seq\\difference_expression\\baseMean_500_fc_1.5\\Filtered_new_' + c1 + '_vs_' + c2 + '.csv' , 'r')
    oldData = csv.reader(f)
    newData = np.loadtxt('D:\\Workspace_New\\data\\RNA\\diff_expression\\basemean_500_fc_1.5\\Filtered_' + c1 + '_' + c2 + '.txt' , skiprows = 1 , 
                         dtype = gene_type)
    n = 0
    for i in islice(oldData , 1 , None):
        gene_name = i[7].lstrip('\'')
        if gene_name in newData['gene_name']:
            n += 1
            
    print n
    
    
        
    