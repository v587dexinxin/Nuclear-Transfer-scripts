# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 17:05:19 2020

@author: han-luo
"""

import numpy as np

data_type = np.dtype({'names':['gene_id' , 'gene_name' , 'FPKM'] , 
                    'formats':['S64' , 'S64' , np.float]})


def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
    

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()           
    
                    
CCS_R1 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\CCS_R1.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
CCS_R2 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\CCS_R2.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
CCS_R3 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\CCS_R3.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
NT5_R1 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\NT5_R1.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
NT5_R2 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\NT5_R2.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
NT5_R3 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\NT5_R3.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
NT5_R4 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\NT5_R4.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
NT6_R1 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\NT6_R1.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
NT6_R2 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\NT6_R2.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
NT6_R3 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\NT6_R3.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
fESC_R1 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\fESC_R1.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
fESC_R2 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\fESC_R2.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))
fESC_R3 = np.loadtxt('H:\\Workspace_New\\data\\RNA\\gene_expression\\fESC_R3.txt' , dtype = data_type , skiprows = 1 , usecols = (0 , 1 , -2))


union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

n = 0
all_gene = []
for i in union_gene:
    gene_name = i['gene_name']
    ccs_1 = CCS_R1[CCS_R1['gene_name'] == gene_name][0]
    ccs_2 = CCS_R2[CCS_R2['gene_name'] == gene_name][0]
    ccs_3 = CCS_R3[CCS_R3['gene_name'] == gene_name][0]
    nt5_1 = NT5_R1[NT5_R1['gene_name'] == gene_name][0]
    nt5_2 = NT5_R2[NT5_R2['gene_name'] == gene_name][0]
    nt5_3 = NT5_R3[NT5_R3['gene_name'] == gene_name][0]
    nt5_4 = NT5_R4[NT5_R4['gene_name'] == gene_name][0]
    nt6_1 = NT6_R1[NT6_R1['gene_name'] == gene_name][0]
    nt6_2 = NT6_R2[NT6_R2['gene_name'] == gene_name][0]
    nt6_3 = NT6_R3[NT6_R3['gene_name'] == gene_name][0]
    fesc_1 = fESC_R1[fESC_R1['gene_name'] == gene_name][0]
    fesc_2 = fESC_R2[fESC_R2['gene_name'] == gene_name][0]
    fesc_3 = fESC_R3[fESC_R3['gene_name'] == gene_name][0]
    all_gene.append((ccs_1['gene_id'] , i['gene_name'] , i['chr'] , i['strand'] , i['start'] , i['end'] , ccs_1['FPKM'], ccs_2['FPKM'], ccs_3['FPKM'], nt5_1['FPKM'], nt5_2['FPKM'], nt5_3['FPKM'] , nt5_4['FPKM'], nt6_1['FPKM'], nt6_2['FPKM'], nt6_3['FPKM'], fesc_1['FPKM'], fesc_2['FPKM'], fesc_3['FPKM']))
    n += 1
    if n // 10000 > 0:
        print i
    
    
Write2fils_nochr('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_genes.txt' , all_gene)    