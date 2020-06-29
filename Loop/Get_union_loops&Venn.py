# -*- coding: utf-8 -*-
"""
Created on Fri Nov 02 22:16:45 2018

@author: xxli
"""

from __future__ import division
import numpy as np
import sys, cPickle
import os
import math
from sklearn.cluster import KMeans



#--------------------------------------------------------function-------------------------------------------------------------
def Load_loop(peak_file):
    """
    """
    loopType = np.dtype({'names':['chr','S','E'],
                     'formats':['S2', np.int, np.int]})

    loop = np.loadtxt(peak_file,dtype = loopType,skiprows = 1, usecols = (0,1,2))
    
    return loop



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

def Get_all_loops(cellname):
    loop_type = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['S8' , np.int , np.int , 'S8']})
    union = []
    for c in cellname:
        f = Get_loops(os.path.join(LoopFolder , c))
        for i in f:
            union.append((i[0] , i[1] , i[2] , c.split("_")[1]))
    union = np.array(union , dtype = loop_type)
    return union

        


def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a



def center_sites(lists):
    sum_x = 0
    sum_y = 0
    for x in lists:
        sum_x += x[1]
        sum_y += x[2]
    n = len(lists)
    return [float(sum_x)/n, float(sum_y)/n]

def distance(sites_1,sites_2):
    dis = math.sqrt((sites_1[0]-sites_2[1])**2 + (sites_1[1]-sites_2[2])**2)
    return dis

def clustering(loop_sites, dis):
    chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19', 'X']
    classes = []
    for i in chrom:
        c_loops = sorted(loop_sites[loop_sites['chr'] == i], key = lambda x:x[1])
        while True :
            c_classs = []
            c_classs.append((c_loops[0][0] , c_loops[0][1] , c_loops[0][2] , c_loops[0][3]))
            c_loops.remove(c_loops[0])
            center = center_sites(c_classs)
            for loop in c_loops:     
                 if distance(center, loop) <= dis:
                      c_classs.append((loop[0] , loop[1] , loop[2] , loop[3]))
                      center = center_sites(c_classs)
                      c_loops.remove(loop)
            classes.append(c_classs)
            if len(c_loops) == 0:
                break
    return classes

def clustering_rescue(clustered_sites , dis , out):
    loop_type = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['S8' , np.int , np.int , 'S8']})
    single = []
    multi = []
    for i in clustered_sites:
        if len(i) > 1:
            multi.append([s for s in i])
        else:
            single.append((i[0][0] , i[0][1] , i[0][2] , i[0][3]))
        
    single = np.array(single , dtype = loop_type)
    for i in range(len(multi)):
        for j in multi[i]:
            chro = j[0]
            s_s = j[1] - dis
            s_e = j[1] + dis
            e_s = j[2] - dis
            e_e = j[2] + dis
            s = single[single['chr'] == chro]
            mask = (s['start'] >= s_s) & (s['start'] <= s_e) & (s['end'] >= e_s) & (s['end'] <= e_e)
            overlap = s[mask]
            if overlap.size != 0:
                single = list(single)
                for k in overlap:
                    clustered_sites[clustered_sites.index(multi[i])].append((k[0] , k[1] , k[2] , k[3]))
                    multi[i].append((k[0] , k[1] , k[2] , k[3]))
                    clustered_sites.remove([(k[0] , k[1] , k[2] , k[3])])
                    single.remove(k)
                single = np.array(single , dtype = loop_type)
            else:
                continue
                

    union_new = [np.array(a , dtype = loop_type) for a in clustered_sites]
    union_sites = []
    for i in union_new:
        chro = i[0][0]
        center = center_sites(i)
        union_sites.append((chro , int(center[0]) , int(center[1])))
    union_sites = Sort(union_sites , 0 , 1)
    o = open(out , 'w')
    for i in union_sites:
        o.writelines('\t'.join([i[0] , str(i[1]) ,str(i[2])]) + '\n')
    o.close()
    return union_new




    
    
LoopFolder = 'D:\\Workspace_New\\data\\HiC\\Loop'
LoopFil1 = 'Clustered_NT6_R1_loops_20K.txt'
LoopFil2 = 'Clustered_NT6_R2_loops_20K.txt'
LoopFil3 = 'Clustered_NT6_loops_20K.txt'
LoopFil4 = 'Clustered_fESC_loops_20K.txt'
R = 20000
out = 'D:\\Workspace_New\\data\\HiC\\Loop\\NT6_R1_R2_20K.txt'
    
initial_distance = 2 * R * math.sqrt(2) + 1000
all_loops = Get_all_loops([LoopFil1 , LoopFil2])
union = clustering(all_loops , initial_distance)
union_new = clustering_rescue(union , 2 * R , out)
	
#-------------------------------------------VennAreaFor4Cell-----------------------------------------------------------------------


n1 = 0 ; n2 = 0 ; n3 = 0 ; n4 = 0 ; n12 = 0 ; n13 = 0 ; n14 = 0 ; n23 = 0 ; n24 = 0 ; n34 = 0 
n123 = 0 ; n124 = 0 ; n134 = 0 ; n234 = 0 ; n1234 = 0
for i in union_new:
    ce = i['cell']
    if 'CCS' in ce:
        n1 += 1
    if 'NT5' in ce:
        n2 += 1
    if 'NT6' in ce:
        n3 += 1
    if 'fESC' in ce:
        n4 += 1
    if ('CCS' in ce) and ('NT5' in ce):
        n12 += 1
    if ('CCS' in ce) and ('NT6' in ce):
        n13 += 1
    if ('CCS' in ce) and ('fESC' in ce):
        n14 += 1
    if ('NT5' in ce) and ('NT6' in ce):
        n23 += 1
    if ('NT5' in ce) and ('fESC' in ce):
        n24 += 1
    if ('NT6' in ce) and ('fESC' in ce):
        n34 += 1
    if ('CCS' in ce) and ('NT5' in ce) and ('NT6' in ce):
        n123 += 1
    if ('CCS' in ce) and ('NT5' in ce) and ('fESC' in ce):
        n124 += 1
    if ('CCS' in ce) and ('NT6' in ce) and ('fESC' in ce):
        n134 += 1
    if ('NT5' in ce) and ('NT6' in ce) and ('fESC' in ce):
        n234 += 1
    if ('CCS' in ce) and ('NT5' in ce) and ('NT6' in ce) and ('fESC' in ce):
        n1234 += 1
    
print 'area1=%d' %n1
print 'area2=%d' %n2
print 'area3=%d' %n3
print 'area4=%d' %n4
print 'n12=%d' %n12
print 'n13=%d' %n13
print 'n14=%d' %n14
print 'n23=%d' %n23
print 'n24=%d' %n24
print 'n34=%d' %n34
print 'n123=%d' %n123
print 'n124=%d' %n124
print 'n134=%d' %n134
print 'n234=%d' %n234
print 'n1234=%d' %n1234   
