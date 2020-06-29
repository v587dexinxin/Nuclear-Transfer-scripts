# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 21:57:44 2018

@author: xxli
"""

import os
import re
import pandas as pd
def merge_csv(filepath):

    filename_list = []
    csv_list = sorted([filepath+'/'+i for i in os.listdir(filepath)])
    for i in csv_list:
        print(i)
    for i in range(len(csv_list)):
        tmp=pd.read_csv(re.sub(r'\\',r'/',csv_list[i]),sep='\t',header=None)
        filename = re.match(r'.*/(.*)\..*',csv_list[i]).group(1)
        tmp.columns = ['gene',filename]
        filename_list.append(filename)
#        tmp = tmp[tmp['gene'] > 'g']
        if i == 0:
            merged_csv = tmp
        else:
            merged_csv=pd.merge(merged_csv,tmp,on='gene',how='outer')

    for i in filename_list:
        merged_csv[i].fillna(0)

    merged_csv.to_csv('{}/merged.csv'.format(filepath),index=False)

merge_csv(r'D:\Workspace_New\data\RNA\reads_count')
