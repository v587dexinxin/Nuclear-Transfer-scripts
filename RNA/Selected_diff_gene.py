# -*- coding: utf-8 -*-
"""
Created on Thu Dec 06 09:49:48 2018

@author: xxli
"""

import numpy as np 
import csv
from itertools import islice


cell = ['CCS_fESC' , 'CCS_NT5' , 'CCS_NT6' , 'fESC_NT5' , 'fESC_NT6' , 'NT5_NT6']



for c in cell:
    c1 = c.split("_")[0] ; c2 = c.split("_")[1]
    
    diffData = open('D:\\Workspace_New\\data\\RNA\\diff_expression\\' + c1 + '_' + c2 + '.txt' , 'r')
    out = open('D:\\Workspace_New\\data\\RNA\\diff_expression\\basemean_500_fc_1.5\\Filtered_' + c1 + '_' + c2 + '.txt' , 'w')
    n = 0
    for i in diffData:
        out.writelines(i)
        break
    for i in islice(diffData , 1 , None):
        a = i.strip().split()
        if 'NA' not in a:
            if (float(a[1]) >= 500) and (float(a[2]) >= 1.5 or float(a[2]) <= -1.5) and float(a[6]) <= 0.05:
                n += 1
                out.writelines(i)
            else:
                continue
        else:
            continue
    print n
    out.close()
    diffData.close()
                
    
    
    