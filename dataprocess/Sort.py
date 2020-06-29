# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 15:26:49 2019

@author: han-luo
"""
import numpy as np

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a
	
def Sort_array(matrix , col):
    '''
    a: list of needed to be sorted
    '''
    index = np.lexsort([matrix[:,col]])
    matrix = matrix[index]
    return matrix