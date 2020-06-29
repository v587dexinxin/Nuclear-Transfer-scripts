# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 11:24:40 2020

@author: han-luo
"""

import numpy as np
import pandas as pd
import csv

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
    


union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_genes.txt')
    
cl = ['CCS_fESC' , 'CCS_NT5' , 'CCS_NT6' , 'fESC_NT5' , 'fESC_NT6' , 'NT5_NT6']
for c in cl:
    out = file('H:\\Workspace_New\\data\\RNA\\diff_expression\\fc_1.5\\Filtered_' + c + '.csv' , 'wb')
    outer = csv.writer(out)
    outer.writerow(['gene_id' , 'baseMean' ,	'log2FoldChange' , 'lfcSE' , 'stat' , 'pvalue' , 'padj' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'fESC_FPKM'])
    df1 = pd.read_csv('H:\\Workspace_New\\data\\RNA\\diff_expression\\' + c + '.csv', index_col=1)
    data = df1[[u'baseMean' , u'log2FoldChange' , u'lfcSE' , u'stat' , u'pvalue' , u'padj']]

    data = data[(data.padj <= 0.05)]
    data = data[(abs(data.log2FoldChange) >= 1.5)]
    Data = []
    for i in range(len(data.index)):
        gene_id = data.index[i]
        genes = union_gene[union_gene['gene_id'] == gene_id]
        if genes.size != 0:
            CCS = (genes[0]['CCS_1'] + genes[0]['CCS_2'] + genes[0]['CCS_3']) / 3
            NT5 = (genes[0]['NT5_1'] + genes[0]['NT5_2'] + genes[0]['NT5_3'] + genes[0]['NT5_4']) / 4
            NT6 = (genes[0]['NT6_1'] + genes[0]['NT6_2'] + genes[0]['NT6_3']) / 3
            fESC = (genes[0]['fESC_1'] + genes[0]['fESC_2'] + genes[0]['fESC_3']) / 3
            if (CCS > 1) or (NT5 > 1) or (NT6 > 1) or (fESC > 1): 
                Data.append((gene_id , str(data.baseMean[i]) , str(data.log2FoldChange[i]) , str(data.lfcSE[i]) , str(data.stat[i]) , str(data.pvalue[i]) , str(data.padj[i]) , genes[0]['gene_name'] , genes[0]['chr'] , genes[0]['strand'] , str(genes[0]['start']) , str(genes[0]['end']) , str(CCS) , str(NT5) , str(NT6) , str(fESC)))
        else:
            pass
    Data = list(set(Data))
    outer.writerows(Data)
    out.close()
