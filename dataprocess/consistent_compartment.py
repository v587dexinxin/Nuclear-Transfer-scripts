# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 09:25:01 2018

@author: xxli
"""


import os
import numpy as np


c = 'fESC_R2'
sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]}) 
PCFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment'
OutFolder = 'D:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new'

PCFil = c + '_Local_Chromosome_Matrix.npz_compartment_200K.txt'
PCSource = os.path.join(PCFolder , PCFil)
PCData = np.loadtxt(PCSource , dtype = sig_type)
OutFil = c + '_compartment_200K.txt'
chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' ,'8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
o = open(os.path.join(OutFolder , OutFil) , 'w')

for g in chrom:
    data = PCData[PCData['chr'] == g]
    if (g == '15') or (g == '16'):
        for i in data:
            o.writelines(i['chr'] + '\t' + str(-i['score']) + '\n')
    else:
        for i in data:
            o.writelines(i['chr'] + '\t' + str(i['score']) + '\n')
             
o.close()
















