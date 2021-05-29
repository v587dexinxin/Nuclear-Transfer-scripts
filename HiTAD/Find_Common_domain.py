# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 20:12:26 2021

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



def Common_domain(tad1 , tad2):
    '''
    '''
    common = []
    for i in tad1:
        chro = i['chr']
        start_s = i['start'] - 40000
        start_e = i['start'] + 40000
        end_s = i['end'] - 40000
        end_e = i['end'] + 40000
        tmp = tad2[tad2['chr'] == chro]
        mask = (tmp['start'] >= start_s) & (tmp['start'] <= start_e) & (tmp['end'] >= end_s) & (tmp['end'] <= end_e)
        overlap = tmp[mask]
        if overlap.size != 0:
            if overlap.size > 1:
                print (i , overlap)
            common.append(i)
    common = np.array(common , dtype = tad1.dtype)
    return common


cell = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']

for c in cell:
    data1 = Load_domain('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\Replicates\\' + c + '_R1_Domain_bottom_40K_respective_stable_600K.txt')
    data2 = Load_domain('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\Replicates\\' + c + '_R2_Domain_bottom_40K_respective_stable_600K.txt')
    common = Common_domain(data1 , data2)
    out = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\Replicates\\' + c + '_R1_R2_common_Domain_bottom_40K_respective_stable_600K.txt' , 'w')
    for i in common:
        out.writelines('\t'.join([i['chr'] , str(i['start']) , str(i['end'])]) + '\n' )
    out.close()
    print (c , len(data1) , len(data2) , len(common))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    