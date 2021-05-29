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
# matplotlib.use('Agg')
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
                     'formats':['U8', np.int, np.int]})

    loop = np.loadtxt(peak_file,dtype = loopType,skiprows = 1, usecols = (0,1,2))
    
    return loop



def Get_loops(LoopFil):
    """
    cell : String ,The name of cells  e.g.'fESC','ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['U8' , np.int , np.int]})
    loops = []
    fileSource = os.path.join(LoopFil)
    Loop = np.loadtxt(fileSource, usecols = (0,1,2) , dtype = loop_type, skiprows = 1)
    for i in Loop:
        if i['end'] - i['start'] >= 300000:
            loops.append((i[0].lstrip('chr') , i[1] , i[2]))
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops



# def Get_all_loops(cellname):
#     loop_type = ({'names':['chr' , 'start' , 'end' , 'cell'],
#                   'formats':['U8' , np.int , np.int , 'U8']})
#     union = []
#     for c in cellname:
#         f = Get_loops(os.path.join(LoopFolder , c))
#         for i in f:
#             union.append((i[0] , i[1] , i[2] , c.split("_")[2]))
#     union = np.array(union , dtype = loop_type)
#     return union


def Get_all_loops(loop_list , cellname):
    loop_type = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['U8' , np.int , np.int , 'U8']})
    union = []
    for i in range(len(loop_list)):
        for j in loop_list[i]:
            union.append((j[0] , j[1] , j[2] , cellname[i]))
    union = np.array(union , dtype = loop_type)
    return union



        
def Get_commpn_loops(loop_1 , loop_2):

    common_1 = []
    common_2 = []
    for i in loop_1:
        chro = i['chr']
        tmp = loop_2[loop_2['chr'] == chro]
        common = 0
        for j in tmp:
            d = distance([i[1] , i[2]] , j)
            if d <= initial_distance:
                common += 1
        if common > 0:
            common_1.append(i)
    for i in loop_2:
        chro = i['chr']
        tmp = loop_1[loop_1['chr'] == chro]
        common = 0
        for j in tmp:
            d = distance([i[1] , i[2]] , j)
            if d <= initial_distance:
                common += 1
        if common > 0:
            common_2.append(i)
        
    common_1 = np.array(common_1 , dtype = loop_1.dtype)
    common_2 = np.array(common_2 , dtype = loop_1.dtype)
    print (len(common_1) , len(common_2))

    return common_1


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
                  'formats':['U8' , np.int , np.int , 'U8']})
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
                  'formats':['U8' , np.int , np.int]})
    common_1 = []
    common_2 = []
    for i in loop_1:
        chro = i['chr']
        tmp = loop_2[loop_2['chr'] == chro]
        common = 0
        for j in tmp:
            d = distance([i[1] , i[2]] , j)
            if d <= 4 * math.sqrt(2) + 0.05 :
                common += 1
        if common > 0:
            common_1.append(i)
    for i in loop_2:
        chro = i['chr']
        tmp = loop_1[loop_1['chr'] == chro]
        common = 0
        for j in tmp:
            d = distance([i[1] , i[2]] , j)
            if d <= 4 * math.sqrt(2) + 0.05:
                common += 1
        if common > 0:
            common_2.append(i)
        
    common_1 = np.array(common_1 , dtype = loop_1.dtype)
    common_2 = np.array(common_2 , dtype = loop_1.dtype)
    
    loop_1 = list(loop_1)
    loop_2 = list(loop_2)
    print (len(loop_1) , len(loop_2), len(common_1) , len(common_2))
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
    print (len(loop_1) , len(loop_2))
    diff_1 = np.array(Sort(loop_1 , 0 , 1 ), dtype = loop_type)
    diff_2 = np.array(Sort(loop_2 , 0 , 1 ), dtype = loop_type)
    return diff_1 , diff_2 , common_1 

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join([i[0].lstrip('chr')] + list(i[1:]))+'\n')
    out.close()            
        
def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
 
    

loop_type = ({'names':['chr' , 'start' , 'end'],
              'formats':['U8' , np.int , np.int]})
R = 20000
initial_distance = 2 * math.sqrt(2) + 0.05


LoopFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new'       
CCSFil = 'Cluster_Selected_CCS_Traditional_Loops_20K_Loops_20K.txt'
NT5Fil = 'Cluster_Selected_NT5_Traditional_Loops_20K_Loops_20K.txt'
NT6Fil = 'Cluster_Selected_NT6_Traditional_Loops_20K_Loops_20K.txt'
F35Fil = 'Cluster_Selected_F35_Traditional_Loops_20K_Loops_20K.txt'
F40Fil = 'Cluster_Selected_F40_Traditional_Loops_20K_Loops_20K.txt'
NTsFil = 'Cluster_Selected_NTs_Traditional_Loops_20K_Loops_20K.txt'
fESCFil = 'Cluster_Selected_fESC_Traditional_Loops_20K_Loops_20K.txt'



out = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\Classify_loops\\union_loops\\NT_fESC_union_loops_20K.txt'


CCS_loop = Get_loops(os.path.join(LoopFolder , CCSFil))
NT5_loop = Get_loops(os.path.join(LoopFolder , NT5Fil))
NT6_loop = Get_loops(os.path.join(LoopFolder , NT6Fil))
F35_loop = Get_loops(os.path.join(LoopFolder , F35Fil))
F40_loop = Get_loops(os.path.join(LoopFolder , F40Fil))
NTs_loop = Get_loops(os.path.join(LoopFolder , NTsFil))
fESC_loop = Get_loops(os.path.join(LoopFolder , fESCFil))


# ##---------------------respective_NT_fESC_union-----------------------

NT5_speci , NT6_speci , NTs = Get_diff_loops(NT5_loop, NT6_loop)
F35_speci , F40_speci , fESC = Get_diff_loops(F35_loop, F40_loop)

all_loops = Get_all_loops([NTs , fESC] , ['NTs' , 'fESC'])
union = clustering(all_loops , initial_distance)
union_new = clustering_rescue(union , 2 * R , out)
NTs_speci , fESC_speci , common_NTs_fESC = Get_diff_loops(NTs , fESC)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\HICCUPS\\union_loops\\NT_fESC_common_loops_20K.txt' , common_NTs_fESC)


CCS_speci , NT_ESC_speci , common = Get_diff_loops(CCS_loop , common_NTs_fESC)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\HICCUPS\\union_loops\\respective_NT_fESC_common\\CCS_specific_loops_20K.txt' , CCS_speci)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\HICCUPS\\union_loops\\respective_NT_fESC_common\\NT_fESC_specific_loops_20K.txt' , NT_ESC_speci)


##-----------------------merge_NT_fESC_union---------------------


all_loops = Get_all_loops([NTs_loop , fESC_loop] , ['NTs' , 'fESC'])
union = clustering(all_loops , initial_distance)
union_new = clustering_rescue(union , 2 * R , out)

NTs_speci , fESC_speci , common_NTs_fESC = Get_diff_loops(NTs_loop , fESC_loop)

Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\Classify_loops\\union_loops\\NT_fESC_common_loops_20K_merge.txt' , common_NTs_fESC)

CCS_speci , NT_ESC_speci , common = Get_diff_loops(CCS_loop , common_NTs_fESC)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\Classify_loops\\union_loops\\CCS_specific_loops_20K_NT_fESC_union.txt' , CCS_speci)
Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Loop_new\\Classify_loops\\union_loops\\NT_fESC_specific_loops_20K_NT_fESC_union.txt' , NT_ESC_speci)



CCS_loop = CCS_loop[CCS_loop['chr'] != 'X']
NTs_loop = NTs_loop[NTs_loop['chr'] != 'X']
fESC_loop = fESC_loop[fESC_loop['chr'] != 'X']








