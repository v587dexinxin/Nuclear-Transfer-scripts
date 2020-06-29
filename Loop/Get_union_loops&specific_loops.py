# -*- coding: utf-8 -*-
"""
Created on Sun Oct 27 18:44:28 2019

@author: han-luo
"""
from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from scipy.special import ndtr
import math
#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')

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
    dis = math.sqrt((sites_1[0] / 20000 - sites_2[1] / 20000)**2 + (sites_1[1] / 20000-sites_2[2] / 20000)**2)
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
    loop_type_1 = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['S8' , np.int , np.int , 'S8']})
    single = []
    multi = []
    for i in clustered_sites:
        if len(i) > 1:
            multi.append([s for s in i])
        else:
            single.append((i[0][0] , i[0][1] , i[0][2] , i[0][3]))
        
    single = np.array(single , dtype = loop_type_1)
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
                single = np.array(single , dtype = loop_type_1)
            else:
                continue
                

    union_new = [np.array(a , dtype = loop_type_1) for a in clustered_sites]
    union_sites = []
    for i in union_new:
        chro = i[0][0]
        center = center_sites(i)
        union_sites.append((chro , int(center[0]) , int(center[1])))
    union_sites = Sort(union_sites , 0 , 1)
    union_sites = np.array(union_sites , dtype = loop_type)
    o = open(out , 'w')
    for i in union_sites:
        o.writelines('\t'.join([i[0] , str(i[1]) ,str(i[2])]) + '\n')
    o.close()
    return union_sites



def Get_diff_loops(loop_1 , loop_2):
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['S8' , np.int , np.int]})
    common_1 = []
    common_2 = []
    for i in loop_1:
        chro = i['chr']
        tmp = loop_2[loop_2['chr'] == chro]
        for j in tmp:
            d = distance([i[1] , i[2]] , j)
            if d <= initial_distance:
                common_1.append(i)
    for i in loop_2:
        chro = i['chr']
        tmp = loop_1[loop_1['chr'] == chro]
        for j in tmp:
            d = distance([i[1] , i[2]] , j)
            if d <= initial_distance:
                common_2.append(i)
    common_1 = np.array(common_1 , dtype = loop_type)
    common_2 = np.array(common_2 , dtype = loop_type)
    loop_1 = list(loop_1)
    loop_2 = list(loop_2)
    print len(loop_1) ; print len(loop_2) ; print len(common_1) ; print len(common_2)
    for i in common_1:
        try:
            loop_1.remove(i)
        except:
            pass
    for j in common_2:
        try:
            loop_2.remove(j)
        except:
            pass
    print len(loop_1) ; print len(loop_2)   
    diff_1 = np.array(Sort(loop_1 , 0 , 1 ), dtype = loop_type)
    diff_2 = np.array(Sort(loop_2 , 0 , 1 ), dtype = loop_type)
    return diff_1 , diff_2 , common_1 , common_2

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
        
        

loop_type = ({'names':['chr' , 'start' , 'end'],
              'formats':['S8' , np.int , np.int]})
R = 20000
initial_distance = 2 * math.sqrt(2) + 0.05


LoopFolder = 'H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05'       
CCSFil = 'Cluster_CCS_loops20K_3.txt'
NT5Fil = 'Cluster_NT5_loops20K_3.txt'
NT6Fil = 'Cluster_NT6_loops20K_3.txt'
fESCFil = 'Cluster_fESC_loops20K_3.txt'

out = 'H:\\Workspace_New\\data\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\NT_fESC_union_loops_20K.txt'


CCS_loop = Get_loops(os.path.join(LoopFolder , CCSFil))
all_loops = Get_all_loops([NT5Fil , NT6Fil , fESCFil])
union = clustering(all_loops , initial_distance)
union_new = clustering_rescue(union , 2 * R , out)

CCS_speci , NT_ESC_speci , common_1 , common_2 = Get_diff_loops(CCS_loop , union_new)
Write2fils_nochr('H:\\Workspace_New\\data\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\CCS_specific_loops_20K.txt' , CCS_speci)
Write2fils_nochr('H:\\Workspace_New\\data\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\NT_fESC_specific_loops_20K.txt' , NT_ESC_speci)
