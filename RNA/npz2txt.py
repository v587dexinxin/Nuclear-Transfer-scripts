# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 10:15:58 2019

@author: han-luo
"""

import numpy as np

matrix = np.load('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\CCS_40K\\Correct_Merged_Reps_Local_Chromosome_Matrix.npz')
matrix = matrix['Matrix'][()]

out = open('H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\CCS_40K\\Correct_Merged_Reps_Local_Chromosome_chr1_Matrix.txt' , 'w')

out.writelines('\t' + '\t'.join(['bin' + str(i + 1) + '|mm10|chr1:' + str(i * 40000) + '-' + str((i + 1) * 40000) for i in range(len(matrix['1']))]) + '\n')
for i in range(len(matrix['1'])):
    out.writelines('bin' + str(i + 1) + '|mm10|chr1:' + str(i * 40000) + '-' + str((i + 1) * 40000) + '\t')
    out.writelines('\t'.join([str(x) for x in matrix['1'][i]]) + '\n')

out.close()






