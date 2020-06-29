# -*- coding: utf-8 -*-
"""
Created on Wed Oct 09 14:47:01 2019

@author: han-luo
"""
import numpy as np
from itertools import islice
def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf



def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene



gene_type = ({'names':['gene_name'],
               'formats':['S64']})

##-----------------------------get_diff_gene_names--------------------------------------------
diff_gene = np.loadtxt('H:\\Workspace_New\\data\\RNA\\diff_expression\\diff_gene_matrix_np.log2(FPKM+1)_classify2.txt' , dtype = gene_type , usecols = (0) , skiprows = 1)
diff_gene = diff_gene['gene_name']  

#------------------------------get_union_gene----------------------------------------------
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')        


CCS_name = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\classify\\CCS_noNT_nofESC_Promoter.bed' , dtype = gene_type , usecols = (-1))
CCS_name = CCS_name['gene_name']
ES_name = np.loadtxt('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\classify\\noCCS_NT_fESC_Promoter.bed' , dtype = gene_type , usecols = (-1))
ES_name = ES_name['gene_name']
CCS_names = []
for i in set(CCS_name):
    i = set(i.split(','))
    for j in i:
        if j in diff_gene:
           CCS_names.append(j)
        else:
            pass


len(CCS_names)

ES_names = []
for i in set(ES_name):
    i = set(i.split(','))
    for j in i:
        if j in diff_gene:
           ES_names.append(j)
        else:
            pass


len(ES_names)

gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')

CCS_id = []
for i in CCS_names:
    gene_id = gtf[gtf['gene_name'] == i]['gene_id'][0]
    CCS_id.append(gene_id)
    
ES_id = []
for i in ES_names:
    gene_id = gtf[gtf['gene_name'] == i]['gene_id'][0]
    ES_id.append(gene_id)
 
w = open('H:\\Workspace_New\\data\\ATAC\\peaks\\idr\\classify\\noCCS_NT_fESC_gene_id.txt' , 'w')
for i in ES_id:
    w.writelines(i + '\n')
w.close()
    