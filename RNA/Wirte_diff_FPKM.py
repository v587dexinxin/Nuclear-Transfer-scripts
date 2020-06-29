# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 15:06:22 2020

@author: han-luo
"""

import csv
import numpy as np
import os 


def get_union_gene(Fil):
    union_type = ({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_1' , 'CCS_2' , 'CCS_3' , 'NT5_1' , 'NT5_2' , 'NT5_3' , 'NT5_4' , 'NT6_1' , 'NT6_2' , 'NT6_3' , 'fESC_1' , 'fESC_2' , 'fESC_3'],
                 'formats':['S64' , 'S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if ((i['CCS_1'] + i['CCS_2'] + i['CCS_3']) / 3 > 1) or ((i['NT5_1'] + i['NT5_2'] + i['NT5_3'] + i['NT5_4']) / 4 > 1) or ((i['NT6_1'] + i['NT6_2'] + i['NT6_3']) / 3 > 1) or ((i['fESC_1'] + i['fESC_2'] + i['fESC_3']) / 3 > 1):
            union.append(i)
    union = np.array(union , dtype = union_type)
    return union_gene
    
    
dataFolder = 'H:\\Workspace_New\\data\\RNA\\diff_expression\\difference_expression'

cl = ['CCS_vs_NT5' , 'NT5_vs_fESC' , 'NT5_vs_NT6']
order = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}

union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_genes.txt')


    

all_genes = {}
for c in cl:
    c1 = c.split('_')[0] ; c2 = c.split('_')[2]
    all_genes[c] = []
    data = csv.reader(open(os.path.join(dataFolder , 'new_' +c + '.csv')))
    if c1 == 'NT5' or c2 == 'NT5':       
        for i in data:
            cell_list = [x.rstrip('_fpkm') for x in [i[11] , i[12] , i[13] , i[14] , i[15] , i[16] , i[17]]]
            d_type = ({'names':['gene_name'] + cell_list,
                       'formats':['S64'] + [np.float for x in range(len(cell_list))]})
            break
        for i in data:
            gene_name = i[7].lstrip("'")
            all_genes[c].append((gene_name , i[11] , i[12] , i[13] , i[14] , i[15] , i[16] , i[17]))  
        all_genes[c] = np.array(all_genes[c] , dtype = d_type)
                
    elif c1 == 'CCS' and c2 == 'NT6':
        for i in data:
            cell_list = ['CCS_1' , 'CCS_2' , 'CCS_3' , 'NT6_1' , 'NT6_2' , 'NT6_3']
            d_type = ({'names':['gene_name'] + cell_list,
                       'formats':['S64'] + [np.float for x in range(len(cell_list))]})
            break
        for i in data:
            gene_name = i[7].lstrip("'")
            all_genes[c].append((gene_name , i[11] , i[12] , i[13] , i[15] , i[16] , i[17]))
        all_genes[c] = np.array(all_genes[c] , dtype = d_type)
        
    else:
        for i in data:
            cell_list = [x.rstrip('_fpkm') for x in [i[11] , i[12] , i[13] , i[14] , i[15] , i[16]]]
            d_type = ({'names':['gene_name'] + cell_list,
                       'formats':['S64'] + [np.float for x in range(len(cell_list))]})
            break
        for i in data:
            gene_name = i[7].lstrip("'")
            all_genes[c].append((gene_name , i[11] , i[12] , i[13] , i[14] , i[15] , i[16]))  
        all_genes[c] = np.array(all_genes[c] , dtype = d_type)
            

for c in cl:
    d_type = all_genes[c].dtype
    cells = list(d_type.names)
    cells.remove('gene_name')
    print c , cells
    for i in all_genes[c]:        
        gene_name = i['gene_name']
        mask = union_gene['gene_name'] == gene_name 
        order = [x for x in range(len(mask)) if mask[x] == True]
        if len(order) == 0:
            print gene_name
        else:
            order = order[0]
            for n in cells:
                union_gene[order][n] = i[n]

            