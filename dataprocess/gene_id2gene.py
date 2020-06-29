# -*- coding: utf-8 -*-
"""
Created on Wed Dec 05 20:53:15 2018

@author: xxli
"""

import numpy as np 
import csv
from itertools import islice

cell = ['CCS_fESC' , 'CCS_NT5' , 'CCS_NT6' , 'fESC_NT5' , 'fESC_NT6' , 'NT5_NT6']

geneFil = 'H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt'
gtfFil = open('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf' , 'r')

gtf_type = ({'names':['gene_id' , 'transcript_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                 'formats':['S64' , 'S64' , 'S64' , 'S8' , 'S8' , np.int , np.int]})

gene_type = ({'names':['gene_name' , 'chr' , 'starnd' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})


gtf = []
geneData = np.loadtxt(geneFil , dtype = gene_type , skiprows = 1)

for i in islice(gtfFil , 5 , None):
    a = i.strip().split()
    if a[2] == 'transcript':
        gene_name = i.strip().split('gene_name')[1].split('\"')[1]
        gene_id = i.strip().split('gene_id')[1].split('\"')[1]
        transcript_id = i.strip().split('transcript_id')[1].split('\"')[1]
        gtf.append((gene_id , transcript_id , gene_name , a[0].lstrip('chr') , a[6] , a[3] , a[4]))

gtf = np.array(gtf , dtype = gtf_type)        
gtfFil.close()


for c in cell:
    cell1 = c.split("_")[0] ; cell2 = c.split("_")[1]
    diffFil = open('D:\\Workspace_New\\data\\RNA\\diff_expression\\' + cell1 + '_' + cell2 + '.csv')
    outFil = open('D:\\Workspace_New\\data\RNA\\diff_expression\\' + cell1 + '_' + cell2 + '.txt' , 'w')
    outFil.writelines('\t'.join(['gene_id' , 'baseMean' , 'log2FoldChange' , 'lfcSE' , 'stat' , 'pvalue' , 'padj' , 'gene_name' , 'chr' , 'strand' , 
                      'start' , 'end' , cell1 + '_FPKM' , cell2 + '_FPKM']) + '\n')
    diffData = csv.reader(diffFil)


    for i in islice(diffData , 1 , None):
        ref_gene = gtf[gtf['gene_id'] == i[1]][0]
        gene_name = ref_gene['gene_name']
        strand = ref_gene['strand']
        gene = geneData[geneData['gene_name'] == gene_name][0]
        outFil.writelines('\t'.join([i[1] , i[2] , i[3] , i[4] , i[5] , i[6] , i[7] , gene['gene_name'] , gene['chr'] , strand , str(gene['start']) ,           
                          str(gene['end']) , str(gene[cell1]) , str(gene[cell2])]) + '\n')
    gtfFil.close()
    outFil.close()