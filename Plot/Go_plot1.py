# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 15:11:58 2018

@author: xxli
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font=matplotlib.font_manager.FontProperties(fname=r"C:\Windows\Fonts\Deng.ttf")

d_type = ({'names':['description'  , 'count'],
          'formats':['S64' , np.int]})
fil1 = 'CCS_up_pathway.txt'
fil2 = 'CCS_down_pathway.txt'

al0 = []
al1 = []

f0 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + fil1 ,'r')
f1 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + fil2 ,'r')
for j in f0:
    j = j.strip().split('\t')
    description = j[1].split('(')[0].strip()
    count = int(j[2])
    al0.append((description , count))
al0 = np.array(al0 , dtype = d_type)
f0.close()
for j in f1:
    j = j.strip().split('\t')
    description = j[1].split('(')[0].strip()
    count = int(j[2])
    al1.append((description , count))
al1 = np.array(al1 , dtype = d_type)
f1.close()

al = []
set0 = list(set(al0['description']).union(set(al1['description'])))

y0 = []
y1 = []
for j in set0:
    mask0 = (al0['description'] == j)
    mask1 = (al1['description'] == j)
    overlap0 = al0[mask0]
    overlap1 = al1[mask1]
    if overlap0.size != 0:
        y0.append(overlap0[0]['count'])
    else:
        y0.append(0)
    if overlap1.size != 0:
        y1.append(overlap1[0]['count'])
    else:
        y1.append(0)
    al = [y0 , y1]
            
    
def Bar_plot(al0 , al1):
    left, bottom, width, height = 0.6, 0.1, 0.3, 0.8
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    x00 = np.arange(0 , (len(set0) + 1) * 3  , 3)
    y00 = [0] + al[0]
    x10 = x00 + 1
    y10 = [0] + al[1]
    ax.barh(x00 , y00 , color = 'tomato' , label = "CCS_up")
    ax.barh(x10 , y10 , color = 'cornflowerblue' , label = "CCS_down")
    yticks = x00
    labels = [''] + list(set0) + ['']
    ax.set_yticks(yticks)
    ax.set_yticklabels(labels,fontsize = 15)
    ax.set_xlabel('Gene Numbers' , fontsize = 20)
    ax.legend()
    
for i in range(len(set0)-1 , -1 , -1):
    if al[0][i] <= 5 and al[1][i] <= 5:
        del set0[i]
        del al[0][i]
        del al[1][i]

#    plt.savefig('C:\\Users\\Administrator\\Desktop\\test.png')