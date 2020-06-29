# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 17:02:39 2020

@author: han-luo
"""

import numpy as np
               
def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
            union.append(i)
    union = np.array(union , dtype = union_type)
    return union
    
    
def get_union_gene_sites(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})
         
                 
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
            gene_site = int((i['start'] + i['end']) / 2)
            union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
    union = np.array(union , dtype = union_type_1)
    return union
    
def get_union_gene_all(Fil):
    union_type = ({'names':['gene_id','gene_name','chr','strand','start','end','CCS_1','CCS_2','CCS_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3','fESC_1','fESC_2','fESC_3'],
                 'formats':['S64' , 'S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if ((i['CCS_1'] + i['CCS_2'] +i['CCS_3']) / 3 > 1) or ((i['NT5_1'] + i['NT5_2'] +i['NT5_3'] + i['NT5_4']) / 4 > 1) or ((i['NT6_1'] + i['NT6_2'] +i['NT6_3']) / 3 > 1) or ((i['fESC_1'] + i['fESC_2'] +i['fESC_3']) / 3 > 1):
            union.append(i)
    union = np.array(union , dtype = union_type)
    return union

def get_raw_genes(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)

    return union_gene
    

    
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

union_gene_site = get_union_gene_sites('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

union_gene_all = get_union_gene_all('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_genes.txt')

