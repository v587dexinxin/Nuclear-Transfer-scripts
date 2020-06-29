# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 22:09:02 2019

@author: han-luo
"""

import numpy as np
import cooler
from scipy import sparse
import sys

DataFile = 'NSC_Maternal_40K.npz'
OutFile = '../NSC_Maternal_40K.npz'
maternal = 'M'



	
def matrix2npz(DataFile , OutFile , maternal):
    Matrix_data = np.load(DataFile)
    tmp = {};Tmp = {}
    for g in Matrix_data.keys():
        if g == 'Y':
            continue
        tmp[maternal + g] = Matrix_data[g]
		
    Tmp['Matrix'] = np.array(tmp)
	
    np.savez_compressed(OutFile, **Tmp)
    
matrix2npz(DataFile , OutFile , maternal)




def matrix2npz(DataFile , OutFile):
    Matrix_data = np.load(DataFile , allow_pickle=True)
    Matrix_data = Matrix_data['Matrix'][()]
    tmp = {}  
    for g in Matrix_data.keys():
        if g == 'Y':
            continue
        tmp[g] = Matrix_data[g]
		
	
    np.savez_compressed(OutFile, **tmp)
    
for i in ['Correct_Merged_Reps' , 'SRR1658716' , 'SRR1658722' , 'SRR1658723' , 'SRR1658724' , 'SRR1658725' , 'SRR1658726']:
    DataFile = i + '_Local_Chromosome_Matrix.npz'
    OutFile = i + '_Local_Chromosome_nonMatrix.npz' 
    matrix2npz(DataFile , OutFile)
    



