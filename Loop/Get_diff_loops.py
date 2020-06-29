# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 12:07:14 2018

@author: xxli
"""

from __future__ import division
import numpy as np
import sys, cPickle
import os
import math
from sklearn.cluster import KMeans

def Get_loops(LoopFil):
    """
    cell : String ,The name of cells  e.g.'fESC','ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['S8' , np.int , np.int]})
    loops = []
    fileSource = os.path.join(LoopFil)
    Loop = np.loadtxt(fileSource, usecols = (0,1,2) , dtype = loop_type, skiprows = 1)
    for i in Loop:
        if i['end'] - i['start'] >= 300000:
            loops.append(i)
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a

def distance(sites_1,sites_2):
    dis = math.sqrt((sites_1[1] / 20000 - sites_2[1] / 20000)**2 + (sites_1[2] / 20000 - sites_2[2] / 20000)**2)
    return dis

def Get_diff_loops(loop_1 , loop_2 , c1 , c2):
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['S8' , np.int , np.int]})
    loop_type_1 = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['S8' , np.int , np.int , 'S8']})
    common_1 = []
    common_2 = []
    for i in loop_1:
        chro = i['chr']
        tmp = loop_2[loop_2['chr'] == chro]
        for j in tmp:
            d = distance(i , j)
            if d <= initial_distance:
                common_1.append(i)
    for i in loop_2:
        chro = i['chr']
        tmp = loop_1[loop_1['chr'] == chro]
        for j in tmp:
            d = distance(i , j)
            if d <= initial_distance:
                common_2.append(i)
    common_1 = np.array(common_1 , dtype = loop_type)
    common_2 = np.array(common_2 , dtype = loop_type)
    loop_1 = list(loop_1)
    loop_2 = list(loop_2)
    print len(loop_1) ; print len(loop_2)
    for i in common_1:
        try:
            loop_1.remove(i)
        except:
            pass
    for j in common_2:
        try:
            loop_2.remove(i)
        except:
            pass
    print len(loop_1) ; print len(loop_2)   
    diff = []
    for i in loop_1:
        diff.append((i['chr'] , i['start'] , i['end'] , c1))
    for i in loop_2:
        diff.append((i['chr'] , i['start'] , i['end'] , c2))
    diff = np.array(Sort(diff , 0 , 1 ), dtype = loop_type_1)
    return diff 

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
        

c1 = 'CCS' ; c2 = 'NT5' ; c3 = 'NT6' ; c4 = 'fESC'
initial_distance = 2 * math.sqrt(2)

loopFil_1 = 'D:\\Workspace_New\\data\\HiC\\Loop\\Loops\\Clustered_loops\\Clustered_' + c1 + '_loops_20K.txt'
loopFil_2 = 'D:\\Workspace_New\\data\\HiC\\Loop\\Loops\\Clustered_loops\\Clustered_' + c2 + '_loops_20K.txt'
loopFil_3 = 'D:\\Workspace_New\\data\\HiC\\Loop\\Loops\\Clustered_loops\\Clustered_' + c3 + '_loops_20K.txt'
loopFil_4 = 'D:\\Workspace_New\\data\\HiC\\Loop\\Loops\\Clustered_loops\\Clustered_' + c4 + '_loops_20K.txt'
unionFil = 'D:\\Workspace_New\\data\\HiC\\Loop\\Loops\\union_sites_loops\\union_sites_20K.txt'


loop_1 = Get_loops(loopFil_1)
loop_2 = Get_loops(loopFil_2)
loop_3 = Get_loops(loopFil_3)
loop_4 = Get_loops(loopFil_4)
union = Get_loops(unionFil)

diff = Get_diff_loops(loop_1 , loop_2 , c1 , c2)
Write2fils_nochr('D:\\Workspace_New\\data\\HiC\\Loop\\Loops\\Diff_loops\\' + c1 + '_' + c2 + '_diff_loops.txt' , diff)

matrix = np.zeros((len(union) ,4))
for i in range(len(union)):
    for j in loop_1:
        d = distance(union[i] , j)
        if d <= initial_distance:
            matrix[i][0] = 1
        else:
            continue
    for j in loop_2:
        d = distance(union[i] , j)
        if d <= initial_distance:
            matrix[i][1] = 1
        else:
            continue
    for j in loop_3:
        d = distance(union[i] , j)
        if d <= initial_distance:
            matrix[i][2] = 1
        else:
            continue
    for j in loop_4:
        d = distance(union[i] , j)
        if d <= initial_distance:
            matrix[i][3] = 1
        else:
            continue
    
out =  open('D:\\Workspace_New\\data\\HiC\\Loop\\Loops\\union_sites_loops\\union.txt' , 'w')
for i in range(len(union)):
    out.writelines('\t'.join([str(x) for x in union[i]]) + '\t' + '\t'.join([str(int(x)) for x in matrix[i]]) + '\n')
out.close()           
