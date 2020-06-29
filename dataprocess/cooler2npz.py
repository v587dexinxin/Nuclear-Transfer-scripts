# -*- coding: utf-8 -*-
"""
Created on Sat May 16 11:37:24 2020

@author: han-luo
"""

import numpy as np
import cooler
from scipy import sparse
import sys


def cooler2npz(DataFile , OutFile , res , maternal):
    Matrix = cooler.Cooler(DataFile + '::' + str(res))
    
    tmp = {}
    for g in chroms:
        if g[0] == maternal: 
            matrix = Matrix.matrix(balance=False).fetch(g)
            tmp[g.lstrip(maternal)] = matrix
		
	
    np.savez_compressed(OutFile, tmp)
    



chromsM = ['M' + i for i in map(str,range(1,20))] + ['MX']
chromsP = ['P' + i for i in map(str,range(1,20))] + ['PX']
chroms = chromsM + chromsP


	

cooler2npz('Merged_Imputated_Haplotype_Multi.cool' , ' , maternal)
