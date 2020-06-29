# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 21:14:03 2020

@author: han-luo
"""

import csv
import os
import numpy as np

union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' ,'NT6' ,'fESC'],
                 'formats':['S64' , 'S8' , np.int ,  np.float , np.float , np.float , np.float]})



diff_folder = 'H:\\Workspace_New\\data\\RNA\\diff_expression\\New_fc_1.5'
boun_folder = 'H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes_new'

cl = ['CCS_noNT_nofESC' , 'CCS_NT_nofESC' , 'noCCS_noNT_fESC' , 'CCS_noNT_fESC',
      'noCCS_NT_fESC' , 'CCS_NT_fESC']

diff_1 = csv.reader(open(os.path.join(diff_folder , 'Filtered_NT5_fESC.csv')))
diff_2 = csv.reader(open(os.path.join(diff_folder , 'Filtered_NT5_fESC.csv')))


diff = []
for i in diff_1:
    gene_name = i[7].lstrip("'")
    diff.append(gene_name)


for i in diff_2:
    gene_name = i[7].lstrip("'")
    diff.append(gene_name)
    
    
diff = list(set(diff))    

cl_diff = {}
for c in cl:
    cl_diff[c] = []
    o = open(os.path.join(boun_folder , c + '_diff_genes.txt') , 'w')
    o.writelines('\t'.join(['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' ,'NT6' ,'fESC']) + '\n')
    data = np.loadtxt(os.path.join(boun_folder , c + '_genes.txt') , skiprows = 1 , dtype = union_type)
    for i in data:
        gene_name = i['gene_name']
        if gene_name in diff:
            cl_diff[c].append(i)
    for i in cl_diff[c]:
        o.writelines('\t'.join(list(np.array(list(i) , dtype = 'S64'))) + '\n')
    o.close()
            
            