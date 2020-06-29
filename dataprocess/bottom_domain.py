# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 16:18:02 2018

@author: xxli
"""

import numpy as np
import os

TADFolder = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\stable_window'

cell = ['CCS' ,'NT5' , 'NT6' ,'fESC']
tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['S4' , np.int , np.int]})

for c in cell:
    OutFil = 'D:\\Workspace_New\\data\\HiC\\HiTAD\\stable_window\\bottom_domain\\' + c + '_bottom_domain_40K.txt'
    o = open(OutFil  ,'w')
    TADFil = c + '_40K_allreps.txt'
    TADSource = os.path.join(TADFolder , TADFil)
    TAData = np.loadtxt(TADSource , dtype = tad_type , usecols = (0 , 1 , 2))
    for i in TAData:
        chro = i['chr']
        start = i['start']
        end = i['end']
        tad = TAData[TAData['chr'] == chro]
        mask = (start <= tad['start']) & (end >= tad['end'])
        overlap = tad[mask]
        if overlap.size == 1:
            o.writelines(i[0] + '\t' + str(i[1]) + '\t' + str(i[2]) + '\n')
        else:
            continue
    o.close()