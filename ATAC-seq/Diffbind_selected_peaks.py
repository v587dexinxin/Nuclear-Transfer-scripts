# -*- coding: utf-8 -*-
"""
Created on Mon Mar 22 14:51:19 2021

@author: xxli
"""
        

import numpy as np

peak_type= np.dtype({'names':['chr' , 'start' , 'end'] , 
                     'formats':['U8' , np.int , np.int]})
diff_type= np.dtype({'names':['chr' , 'start' , 'end' , 'fc' , 'p_value' , 'q_value'] , 
                     'formats':['U8' , np.int , np.int , np.float , np.float , np.float]})


data1 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\Cluster3.bed' , dtype = peak_type)
data2 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\DiffBind\\NTs_vs_fESC_edgeR.txt' , dtype = diff_type , skiprows=1 , usecols=(0 , 1 , 2 , -3 , -2 , -1))

data2 = data2[data2['p_value'] <= 0.05]


repro = [] ; partial = [] ; over = []
for i in data1:
    tmp2 = data2[data2['chr'] == i['chr']]
    mask = (i['start'] <= tmp2['end']) & (i['end'] >= tmp2['start'])
    overlap = tmp2[mask]
    if overlap.size != 0:        
        if overlap[0]['fc'] < 0:
            partial.append(i)
        elif overlap[0]['fc'] > 0:
            over.append(i)
        else:
            repro.append(i)
    else:
        repro.append(i)
            
out1 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\Cluster3_Partial.bed' , 'w')
for i in partial:
    out1.writelines('\t'.join([str(x) for x in i]) + '\n')
out1.close()            
        
out2 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\Cluster3_Over.bed' , 'w')
for i in over:
    out2.writelines('\t'.join([str(x) for x in i]) + '\n')
out2.close()     

out3 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\Cluster3_Repro.bed' , 'w')       
for i in repro:
    out3.writelines('\t'.join([str(x) for x in i]) + '\n')       
out3.close()
        
print (len(repro) , len(partial) , len(over))
        


##-------------------------Cluster2--------------------------------------------


data1 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\DiffBind\\Cluster2.bed' , dtype = peak_type)
data2 = np.loadtxt('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\DiffBind_1\\Cluster2_NTs_vs_fESC_edgeR.txt' , dtype = diff_type , skiprows=1 , usecols=(0 , 1 , 2 , -3 , -2 , -1))

data2 = data2[data2['p_value'] <= 0.05]


repro = [] ; partial = [] ; over = []
n = 0
for i in data1:
    tmp2 = data2[data2['chr'] == i['chr']]
    mask = (i['start'] <= tmp2['end']) & (i['end'] >= tmp2['start'])
    overlap = tmp2[mask]
    if overlap.size != 0:    
        n += 1
        if overlap[0]['fc'] < -0.5:
            over.append(i)
        elif overlap[0]['fc'] > 0.5:
            partial.append(i)
    else:
        repro.append(i)
            
out1 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\DiffBind_1\\Cluster2_partial.bed' , 'w')
for i in partial:
    out1.writelines('\t'.join([str(x) for x in i]) + '\n')
out1.close()            
        
out2 = open('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\ATAC\\ATAC_new\\peaks\\idr_new\\classify_peaks\\DiffBind_1\\Cluster2_over.bed' , 'w')
for i in over:
    out2.writelines('\t'.join([str(x) for x in i]) + '\n')
out2.close()            
 




