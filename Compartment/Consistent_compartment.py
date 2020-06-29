# -*- coding: utf-8 -*-
"""
Created on Thu Nov 07 10:56:52 2019

@author: han-luo
"""

import numpy as np


sig_type = np.dtype({'names':['chr','score'],
                      'formats':['S4',np.float]}) 
                      
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = sig_type )
MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\MEF_PC_200K\\MEF_PC_200K_Compartment_200K.txt' , dtype = sig_type)

fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = sig_type )
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\E14_PC_200K\\E14_PC_200K_Compartment_200K.txt' , dtype = sig_type)
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\IPS_P3_PC_200K\\IPS_P3_PC_200K_Compartment_200K.txt' , dtype = sig_type)

def consitent_compartment(PC1 , PC2 , out):
    '''
    PC1: original compartment
    PC2: the compartment that may not be consitent with PC1
    '''
    chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
    o = open(out , 'w')
    for g in chrom:
        tmp1 = PC1[PC1['chr'] == g]
        tmp2 = PC2[PC2['chr'] == g]
        n = 0
        for i in range(len(tmp1)):
            if tmp1[i]['score'] * tmp2[i]['score'] < 0:
                n += 1
        print n / len(tmp1)
        for i in tmp2:
            if n / len(tmp1) > 0.5:
                o.writelines('\t'.join([i['chr'] , str(-i['score'])]) + '\n')
            else:
                o.writelines('\t'.join([i['chr'] , str(i['score'])]) + '\n')
    o.close()

consitent_compartment(CCS , MEF , 'H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_Compartment_200K.txt')            
consitent_compartment(fESC , E14 , 'H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_Compartment_200K.txt')
consitent_compartment(fESC , IPS_P3 , 'H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_Compartment_200K.txt')

    

