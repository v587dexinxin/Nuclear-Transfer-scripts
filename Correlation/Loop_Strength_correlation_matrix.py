# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 20:20:44 2018

@author: xxli
"""

import numpy as np


union_type = np.dtype({'names':['chr','start','end','CCS','NT5','NT6','fESC'],
                       'formats':['S2' , np.int , np.int , np.float , np.float , np.float , np.float]})
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
ce = ['CCS_NT5' ,'CCS_NT6' , 'CCS_fESC' ,'NT5_NT6' , 'NT5_fESC' , 'NT6_fESC']

cor_matrix = np.zeros((4,4))
for i in ce:
    f = 'D:\\Workspace_New\\data\\HiC\\Loop\\' + i + '_ave.txt'
    loop = np.loadtxt(f , dtype = union_type , skiprows = 1)
    c1 = i.split("_")[0]
    c2 = i.split("_")[1]
    x = []
    y = []
    for j in range(len(loop[c1])):
        if (loop[c1][j] < 50) and (loop[c2][j] < 50):
            x.append(loop[c1][j])
            y.append(loop[c2][j])
    print np.percentile(x , 99.99) , np.percentile(y , 99.99)
    
    cor_matrix[cells[c1]][cells[c2]] = round(np.corrcoef(x , y)[0][1] , 3)
    cor_matrix[cells[c2]][cells[c1]] = round(np.corrcoef(x , y)[0][1] , 3)

cor_matrix[0,0] = 1
cor_matrix[1,1] = 1
cor_matrix[2,2] = 1
cor_matrix[3,3] = 1
 
loop_all = {'CCS':[] , 'NT5':[] , 'NT6':[] , 'fESC':[]}

for i in loop:
    loop_all['CCS'].append(i['CCS'])
    loop_all['NT5'].append(i['NT5'])
    loop_all['NT6'].append(i['NT6'])
    loop_all['fESC'].append(i['fESC'])
    
    
for k , v in loop_all.items():
    v = np.array(v)
    loop_all[k] = v
    
cor_matrix = np.zeros((4,4))

for i in cell:
    for j in cell:
        cor_matrix[cells[i]][cells[j]] = round(np.corrcoef(loop_all[i] , loop_all[j])[0][1] , 3)
        
o = open('D:\\Workspace_New\\data\\HiC\\correlation\\union_in_pairs_loop_strength_correlation_ave.txt' ,  'w')
o.writelines('\t'+'\t'.join(['CCS','NT5','NT6','fESC']) + '\n')
for c in cell:
    o.writelines(c + '\t' + '\t'.join([str(x) for x in cor_matrix[cells[c]]]) + '\n')
o.close()
