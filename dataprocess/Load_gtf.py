# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 16:49:23 2018

@author: xxli
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

gtf = Load_gtf('G:\\data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')

interval = []
for i in ['NANOG' , 'ESRRB' , 'POU5F1' , 'SOX2' , 'KLF4' , 'FOSL1' , 'FOSL2' , 'RUNX1']:
    mask = (gtf['gene_name'] == i)
    overlap = gtf[mask]
    interval.append((overlap[0]['chr'] , overlap[0]['start'] , overlap[0]['end'] , i))
    
o = open('D:\\Workspace_New\\human\\Selected_intervals_hg19.txt' , 'w')
o.writelines('\t'.join(['gene_name' , 'chr' , 'start' , 'end']) + '\n')
for i in interval:
    o.writelines('\t'.join([i[3] , 'chr' + i[0] , str((i[1] // 40000) * 40000 - 1000000) , str((i[1] // 40000) * 40000 + 1000000) ]) + '\n')
o.close()