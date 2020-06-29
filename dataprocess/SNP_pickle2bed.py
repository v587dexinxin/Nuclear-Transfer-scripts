# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 10:14:29 2019

@author: han-luo
"""

import pickle

SnpFil = pickle.load(open('C:\\Users\\lenovo\\Desktop\\PICTURE_New\\SNPS\\ICM_Snps.pickle' , 'rb'))
OutFil = open('C:\\Users\\lenovo\\Desktop\\PICTURE_New\\SNPS\\ICM_Snps.bed' , 'w')


OutFil.writelines('\t'.join(['SNP' , 'Chromosome' , 'Position']) + '\n')
n = 0
for g in [str(x) for x in range(1,20)] + ['X' , 'Y']:
    for i in SnpFil[g]['pos']:
        OutFil.writelines('\t'.join(['snp_' + str(n) , g , str(i)]) + '\n')
        n += 1
OutFil.close()