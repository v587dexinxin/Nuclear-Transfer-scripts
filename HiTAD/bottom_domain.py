# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 11:12:55 2020

@author: xxli
"""


import numpy as np


def Load_domain(tad_fil):
    '''
    '''
    data_type = np.dtype({'names':['chr' , 'start' , 'end'],
                          'formats':['U8' , np.int , np.int]})
    
    data = np.loadtxt(tad_fil , usecols=(0 , 1 , 2) , dtype = data_type)
    
    return data


def Get_bottom_domain(tad):
    '''
    '''
    tad_bottom = []
    for g in chrom:
        tmp = tad[tad['chr'] == g]
        tad_bottom.append(tmp[0])
        for i in range(1 , len(tmp)):
            start = tmp[i]['start']
            if start >= tad_bottom[-1][-1]:
                tad_bottom.append(tmp[i])
            else:
                pass
            
    tad_bottom = np.array(tad_bottom , dtype = tad.dtype)
    
    
    return tad_bottom
                
                
def Write2fils(filname , tad):    
    '''
    '''
    with open(filname,'w') as out:
        for i in tad:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
                    
                
                           
                

cell = ['NTs' , 'fESC']
chrom = [str(x) for x in range(1 , 20)] + ['X']
for c in cell:
    tad_fil = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\' + c + '_Domain_40K_respective_stable_400K.txt'
    tad = Load_domain(tad_fil)
    tad_bottom = Get_bottom_domain(tad)
    out_fil = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\' + c + '_Domain_bottom_40K_respective_stable_400K.txt'
    Write2fils(out_fil , tad_bottom)
    
    print (len(tad) , len(tad_bottom))
    









