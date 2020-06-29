# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 09:45:44 2018

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

al0 = {'MF':[],'BP':[],'CC':[]}
al1 = {'MF':[],'BP':[],'CC':[]}
al2 = {'MF':[],'BP':[],'CC':[]}
for i in ['MF' , 'BP' , 'CC']:
    f0 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_0 + '_' + i + '.txt' ,'r')
    f1 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_1 + '_' + i + '.txt' ,'r')
    f2 = open('D:/Workspace_New/data/RNA/diff_expression/GO/classify_new/' + classify_2 + '_' + i + '.txt' ,'r')
    for j in f0:
        j = j.strip().split('\t')
        description = j[1].split('(')[0].strip()
        count = int(j[2])
        al0[i].append((description , count))
    al0[i] = np.array(al0[i] , dtype = d_type)
    f0.close()
    for j in f1:
        j = j.strip().split('\t')
        description = j[1].split('(')[0].strip()
        count = int(j[2])
        al1[i].append((description , count))
    al1[i] = np.array(al1[i] , dtype = d_type)
    f1.close()
    for j in f2:
        j = j.strip().split('\t')
        description = j[1].split('(')[0].strip()
        count = int(j[2])
        al2[i].append((description , count))
    al2[i] = np.array(al2[i] , dtype = d_type)
    f2.close()

al = {'MF':[] , 'BP':[] , 'CC':[]}
set0 = set(al0['MF']['description']).union(set(al1['MF']['description'])).union(set(al2['MF']['description']))
set1 = set(al0['BP']['description']).union(set(al1['BP']['description'])).union(set(al2['BP']['description']))
set2 = set(al0['CC']['description']).union(set(al1['CC']['description'])).union(set(al2['CC']['description']))
sets = {'MF':set0 , 'BP':set1 , 'CC':set2}

for i in ['MF' , 'BP' , 'CC']:
    y0 = []
    y1 = []
    y2 = []
    for j in sets[i]:
        mask0 = (al0[i]['description'] == j)
        mask1 = (al1[i]['description'] == j)
        mask2 = (al2[i]['description'] == j)
        overlap0 = al0[i][mask0]
        overlap1 = al1[i][mask1]
        overlap2 = al2[i][mask2]
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
    al[i] = [y0 , y1 , y2]
            
    
def Bar_plot(al0 , al1 , al2):
    left, bottom, width, height = 0.5, 0.1, 0.4, 0.8
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    x00 = np.arange(0 , (len(set0) + 1) * 3  , 3)
    x01 = np.arange(0 , len(set1) * 3 , 3) +  x00[-1] + 10 
    x02 = np.arange(0 , (len(set2) + 1) * 3 , 3) +  x01[-1] + 10 
    y00 = [0] + al['MF'][0]
    y01 = al['BP'][0]
    y02 = al['CC'][0] + [0]
    x10 = x00 + 1
    x11 = x01 + 1
    x12 = x02 + 1
    y10 = [0] + al['MF'][1]
    y11 = al['BP'][1]
    y12 = al['CC'][1] + [0]
    x20 = x00 + 2
    x21 = x01 + 2
    x22 = x02 + 2
    y20 = [0] + al['MF'][2]
    y21 = al['BP'][2]
    y22 = al['CC'][2] + [0]
    ax.barh(x00 , y00 , color = 'limegreen' , label = "classify3")
    ax.barh(x01 , y01 , color = 'limegreen')
    ax.barh(x02 , y02 , color = 'limegreen')
    ax.barh(x10 , y10 , color = 'orangered' , label = "classify4")
    ax.barh(x11 , y11 , color = 'orangered')
    ax.barh(x12 , y12 , color = 'orangered')
    ax.barh(x20 , y20 , color = 'orchid' , label = "classify5")
    ax.barh(x21 , y21 , color = 'orchid')
    ax.barh(x22 , y22 , color = 'orchid')
    yticks = np.array(list(x10) + list(x11) + list(x12))
    labels = [''] + list(set0) + list(set1) + list(set2) + ['']
    ax.set_yticks(yticks)
    ax.set_yticklabels(labels,fontsize = 15)
    ax.set_xlabel('Gene Numbers' , fontsize = 20)
    ax.legend()
#    plt.savefig('C:\\Users\\Administrator\\Desktop\\test.png')