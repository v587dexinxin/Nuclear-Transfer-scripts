# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 21:27:54 2018

@author: xxli
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font=matplotlib.font_manager.FontProperties(fname=r"C:\Windows\Fonts\Deng.ttf")

d_type = ({'names':['description'  , 'count'],
          'formats':['S64' , np.int]})
classify_0 = 'c3'
classify_1 = 'c4'
classify_2 = 'c5'

al0 = []
al1 = []
al2 = []

f0 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_0 + '_pathway.txt' ,'r')
f1 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_1 + '_pathway.txt' ,'r')
f2 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_2 + '_pathway.txt' ,'r')
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
for j in f2:
    j = j.strip().split('\t')
    description = j[1].split('(')[0].strip()
    count = int(j[2])
    al2.append((description , count))
al2 = np.array(al2 , dtype = d_type)
f2.close()

al = []
set0 = list(set(al0['description']).union(set(al1['description'])).union(set(al2['description'])))
sets = {'pathway':set0}

for i in sets:
    y0 = []
    y1 = []
    y2 = []
    for j in sets[i]:
        mask0 = (al0['description'] == j)
        mask1 = (al1['description'] == j)
        mask2 = (al2['description'] == j)
        overlap0 = al0[mask0]
        overlap1 = al1[mask1]
        overlap2 = al2[mask2]
        if overlap0.size != 0:
            y0.append(overlap0[0]['count'])
        else:
            y0.append(0)
        if overlap1.size != 0:
            y1.append(overlap1[0]['count'])
        else:
            y1.append(0)
        if overlap2.size != 0:
            y2.append(overlap2[0]['count'])
        else:
            y2.append(0)
    al = [y0 , y1 , y2]
            
    
def Bar_plot(al0 , al1 , al2):
    left, bottom, width, height = 0.6, 0.1, 0.3, 0.8
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    x00 = np.arange(0 , (len(set0) + 1) * 3  , 3)
    y00 = [0] + al[0]
    x10 = x00 + 1
    y10 = [0] + al[1]
    x20 = x00 + 2
    y20 = [0] + al[2]
    ax.barh(x00 , y00 , color = 'cornflowerblue' , label = "classify3")
    ax.barh(x10 , y10 , color = 'tomato' , label = "classify4")
    ax.barh(x20 , y20 , color = 'greenyellow' , label = "classify5")
    yticks = x10
    labels = [''] + list(set0) + ['']
    ax.set_yticks(yticks)
    ax.set_yticklabels(labels,fontsize = 13)
    ax.set_xlabel('Gene Numbers' , fontsize = 20)
    ax.legend()
    
for i in range(len(set0)-1 , -1 , -1):
    if al[0][i] <= 2 and al[1][i] <= 2 and al[2][i] <= 2:
        del set0[i]
        del al[0][i]
        del al[1][i]
        del al[2][i]
#    plt.savefig('C:\\Users\\Administrator\\Desktop\\test.png')