# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 15:48:19 2018

@author: xxli
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

font=matplotlib.font_manager.FontProperties(fname=r"C:\Windows\Fonts\Deng.ttf")

d_type = ({'names':['description'  , 'count'],
          'formats':['S64' , np.int]})

classify_1 = 'c1'
classify_2 = 'c2'

classify_3 = 'c4'
classify_4 = 'c5'

al1 = []
al2 = []
al3 = []
al4 = []

f1 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_1 + '_pathway.txt' ,'r')
f2 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_2 + '_pathway.txt' ,'r')
f3 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_3 + '_pathway.txt' ,'r')
f4 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_4 + '_pathway.txt' ,'r')

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
for j in f3:
    j = j.strip().split('\t')
    description = j[1].split('(')[0].strip()
    count = int(j[2])
    al3.append((description , count))
al3 = np.array(al3 , dtype = d_type)
f3.close()
for j in f4:
    j = j.strip().split('\t')
    description = j[1].split('(')[0].strip()
    count = int(j[2])
    al4.append((description , count))
al4 = np.array(al4 , dtype = d_type)
f4.close()

al = []
set0 = list(set(al1['description']).union(set(al2['description'])).union(set(al3['description'])).union(set(al4['description'])))
sets = {'pathway':set0}

y1 = []
y2 = []
y3 = []
y4 = []
for j in set0:
    mask1 = (al1['description'] == j)
    mask2 = (al2['description'] == j)
    mask3 = (al3['description'] == j)
    mask4 = (al4['description'] == j)
    overlap1 = al1[mask1]
    overlap2 = al2[mask2]
    overlap3 = al3[mask3]
    overlap4 = al4[mask4]
    if overlap1.size != 0:
        y1.append(overlap1[0]['count'])
    else:
        y1.append(0)
    if overlap2.size != 0:
        y2.append(overlap2[0]['count'])
    else:
        y2.append(0)
    if overlap3.size != 0:
        y3.append(overlap3[0]['count'])
    else:
        y3.append(0)
    if overlap4.size != 0:
        y4.append(overlap4[0]['count'])
    else:
        y4.append(0)
    al = [y1 , y2 , y3 , y4]
            
    
def Bar_plot(al0 , al1 , al2):
    left, bottom, width, height = 0.6, 0.1, 0.3, 0.8
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    x00 = np.arange(0 , (len(set0) + 1) * 5  , 5)
    y00 = [0] + al[0]
    x10 = x00 + 1
    y10 = [0] + al[1]
    x20 = x00 + 2
    y20 = [0] + al[2]
    x30 = x00 + 3
    y30 = [0] + al[3]
    ax.barh(x00 , y00 , color = 'crimson' , label = "CCS_up_partical")
    ax.barh(x10 , y10 , color = 'limegreen' , label = "CCS_up_over")
    ax.barh(x20 , y20 , color = 'royalblue' , label = "CCS_down_partical")
    ax.barh(x30 , y30 , color = 'yellow' , label = "CCS_down_over")
    yticks = x10
    labels = [''] + list(set0) + ['']
    ax.set_yticks(yticks)
    ax.set_yticklabels(labels,fontsize = 13)
    ax.set_xlabel('Gene Numbers' , fontsize = 20 )
    ax.legend(loc = 'upper right')
    
for i in range(len(set0)-1 , -1 , -1):
    if al[0][i] <= 1 and al[1][i] <= 1 and al[2][i] <= 1 and al[3][i] <= 1:
        del set0[i]
        del al[0][i]
        del al[1][i]
        del al[2][i]
        del al[3][i]
#    plt.savefig('C:\\Users\\Administrator\\Desktop\\test.png')