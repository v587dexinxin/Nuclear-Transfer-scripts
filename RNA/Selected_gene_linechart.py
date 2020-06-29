# -*- coding: utf-8 -*-
"""
Created on Fri Apr 03 17:14:37 2020

@author: han-luo
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from palettable.colorbrewer.qualitative import  Paired_10
color = Paired_10.hex_colors
#color = color[:2]

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

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
  
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')


#gene_list = ['Magea2', 'Magea3', 'Magea5', 'Magea6', 'Magea8', 'Rhox1', 'Rhox13' , 'Xlr']
#gene_list = ['Dnmt3a', 'Dnmt3l']
gene_list = ['Trp53' , 'Pdk1' , 'Rad51' , 'Sall4' , 'Rad50' , 'Brca1' , 'Parp1']


fig = plt.figure(figsize = (12, 12))
for i , gene in enumerate(gene_list):
    tmp = union_gene[union_gene['gene_name'] == gene][0]
    tmp_ccs = tmp['CCS']
    tmp_nt5 = tmp['NT5']
    tmp_nt6 = tmp['NT6']
    tmp_fesc = tmp['fESC']
       
    plt.plot(range(4),[np.log2(tmp_ccs+1), np.log2(tmp_nt5+1), np.log2(tmp_nt6+1), np.log2(tmp_fesc+1)], c = color[i], marker = 'o', label = gene)
    plt.xticks(range(4), ['CCS','NT5','NT6','fESC'])
    plt.legend(loc = 'upper left')
    plt.ylabel('Log2(fpkm+1)')


run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S3_figs\\DNA_Damage_repair_genes_linechart.pdf')
