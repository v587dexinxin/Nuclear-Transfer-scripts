# -*- coding: utf-8 -*-
"""
Created on Fri Jan 05 14:38:44 2018

@author: xxli
"""


import numpy as np
from scipy import sparse


chrom=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','X']


    
def matrix2sparse(filename , OutFile):
    data = np.load(filename , allow_pickle=True)
    data_new = {}
    for g in chrom:
        
        matrix = sparse.triu(sparse.coo_matrix(data[g]))
        data_new[g] = np.array(zip(matrix.row,matrix.col,matrix.data) , dtype=[ ('bin1', '<i4'), ('bin2', '<i4'), ('IF', '<f8')])
    np.savez_compressed(OutFile, **data_new)

    
matrix2sparse('Correct_Merged_Reps_Local_Chromosome_nonMatrix.npz' , 'Correct_Merged_Reps_Local_Chromosome_Matrix_sparse.npz')