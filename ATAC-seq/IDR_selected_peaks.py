# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:04:57 2021

@author: xxli
"""

import numpy as np


def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a



peak_type = np.dtype({'names':['chr1' , 'start1' , 'end1' , 'name1' , 'score1' , 'strand1' , 'signal1' , 'p_value1' , 'q_value1' , 'peak1' , 'chr2' , 'start2' , 'end2' , 'name2' , 'score2' , 'strand2' , 'signal2' , 'p_value2' , 'q_value2' , 'peak2'],
                      'formats':['U8' , np.int , np.int , 'U64' , np.float , 'U8' , np.float , np.float , np.float , np.float , 'U8' , np.int , np.int , 'U64' , np.float , 'U8' , np.float , np.float , np.float , np.float]})

selected_type = np.dtype({'names':['score1' , 'score2'],
                          'formats':[np.float , np.float]})



cells = ['F40']
for c in cells:
    peaks = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\' + c + '_R1_R2_common_peaks.narrowPeak' , dtype = peak_type)
    idr_selected = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\' + c + '_R1_R2_selected_1.txt' , dtype = selected_type)
    out = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\' + c + '_IDR_Selected_1.narrowPeak' , 'w')
    out1 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\' + c + '_IDR_Selected_1.bed' , 'w')
    
    selected_peaks = []
    
    for i in idr_selected:
        overlap = peaks[(peaks['p_value1'] == i['score1']) & (peaks['p_value2'] == i['score2'])][0]
        selected_peaks.append(overlap)
    selected_peaks = Sort(selected_peaks , 0 , 1)
    print (len(selected_peaks))
    for j in selected_peaks:
        out.writelines('\t'.join([str(x) for  x in j]) + '\n')
        out1.writelines('\t'.join([str(x) for x in list(j)[:3]]) + '\n')
        
    out.close()
    out1.close()
    
        
    
    

cells = ['CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']
for c in cells:
    peaks = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\' + c + '_IDR_Selected.narrowPeak' , dtype = peak_type)
    out = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\' + c + '_IDR_Selected_1.narrowPeak' , 'w')

    for i in peaks:
        out.writelines('\t'.join(str(x) for x in [i['chr1'] , i['start1'] , i['end1'] , i['name1'] , i['score1'] , i['strand1'] , i['signal1'] , i['p_value1'] , i['q_value1'] , i['peak1']]) + '\n')
    out.close()
    














    