# -*- coding: utf-8 -*-
"""
Created on Wed Mar 04 15:07:02 2020

@author: han-luo
"""

import os
import os.path
import numpy as np
import csv

gene_type = np.dtype({'names': ['gene_name'],
                      'formats':['S64']})

DataFolder = 'H:\\Workspace_New\\data\\RNA\\diff_expression\\New_fc_1.5'
file_list = ['Filtered_NT5_fESC.csv' , 'Filtered_NT6_fESC.csv'] 

diff_gene = []

for f in file_list:
    gene = csv.reader(open(os.path.join(DataFolder , f) , 'r'))
    for i in gene:
        diff_gene.append(i[7].lstrip("'"))

diff_gene = set(diff_gene)



cluster_gene = {}

for i in range(1,9):
    cluster_gene[i] = []
    gene = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster' + str(i) + '_geneName.txt' , dtype = gene_type)
    for j in gene:
        if j[0] in diff_gene:
            cluster_gene[i].append(j[0])
        else:
            pass




for k , v in cluster_gene.items():
    print k , len(v)
    o = open('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster' + str(k) + '_diff_geneName_1.txt' , 'w')
    for i in v:
        o.writelines(i + '\n')
    o.close()

        