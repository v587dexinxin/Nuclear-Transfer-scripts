# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 18:13:30 2019

@author: han-luo
"""

import os

data_path = 'F:\\SNP_Sampling\\SNP_900W'
Filenames = os.listdir(data_path)

for filename in Filenames:
    data = open(os.path.join(data_path, filename) , 'r')
    out = open(os.path.join(data_path, filename.split('.')[0] + '.bed') , 'w')
    out.writelines('\t'.join(['SNP' , 'Chromosome' , 'Position']) + '\n')
    n = 0
    for i in data:
        i = i.strip().split()
        out.writelines('\t'.join(['snp_' + str(n) , i[0] , i[1]]) + '\n')
        n += 1
    data.close()
    out.close()


