# -*- coding: utf-8 -*-
"""
Created on Tue Dec 04 15:48:32 2018

@author: xxli
"""


import numpy as np

cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']

tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})

for c in cell:
    TAData = np.loadtxt('D:/Workspace_New/data/HiC/HiTAD/respective/stable_680K/' + c + '_40K_allreps.txt' , dtype = tad_type , usecols = (0 , 1 , 2))
    out = open('D:/Workspace_New/data/HiC/HiTAD/respective/stable_680K/bottom_domain/' + c + '_40K_allreps.txt' , 'w')
    new = []
    for i in TAData:
        chro = i['chr']
        start = i['start']
        end = i['end']
        tad = TAData[TAData['chr'] == chro]
        mask = (start <= tad['start']) & (end >= tad['end'])
        overlap = tad[mask]
        if overlap.size <= 1:
            new.append(i)
        else:
            continue
    new.sort(key = lambda x:(x[0] , x[1]))
    new_1 = []
    for i in range(len(new) - 1):
        if (new[i + 1][1] - new[i][2] >= 3 * 40000) and (new[i + 1][1] - new[i][2] <= 5 * 40000):
            new_1.append((new[i][0] , new[i][2] , new[i + 1][1]))
    new = new + new_1
    new.sort(key = lambda x:(x[0] , x[1]))
    
    for i in new:
        out.writelines('\t'.join([i[0] , str(i[1]) , str(i[2])]) + '\n')
    out.close()
            