# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 10:27:28 2018

@author: xxli

Useage: python union_loop_strength_npz.py  Resolution LoopFil loop_strength_point.txt loop_strength_ave.txt

"""



from __future__ import division
import numpy as np
import sys, cPickle
import os
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



def Get_Loop_Strength(inter,start,end,Len):
    """
    """
    maskA = (inter['bin1'] == start) & (inter['bin2'] == end)
    maskB = (inter['bin1'] >= start - Len) & (inter['bin1'] <= start + Len) & \
            (inter['bin2'] >= end - Len) & (inter['bin2'] <= end + Len)
#    maskC = (inter['bin1'] >= start - 4 * Len) & (inter['bin1'] <= start - 2 * Len) & \
#            (inter['bin2'] >= end + 2 * Len) & (inter['bin2'] <= end + 4 * Len)
    maskD = (inter['bin1'] >= start + 2 * Len) & (inter['bin1'] <= start + 4 * Len) & \
            (inter['bin2'] >= end - 4 * Len) & (inter['bin2'] <= end - 2 * Len)

    A = inter[maskA]
    B = inter[maskB]
#    C = inter[maskC]
    D = inter[maskD]
    n_B = len(np.nonzero(B)[0])
    n_D = len(np.nonzero(D)[0])
    if A.size != 0:
        A_S = float(A['IF'])
    else:
        A_S = 0
    if B.size != 0:
        B_S = B['IF'].sum() / n_B
    else:
        B_S = 0
#        C_S = C['IF'].sum() / ((Len + 2) **2)
    if D.size != 0:
        D_S = D['IF'].sum() / n_D
    else:
        D_S = 1
    return A_S /D_S ,B_S /D_S


def read_npz(chro,HiCData):
    """
    HiCData in every chro
    """
    data = {}
    inter = HiCData[chro]
    max_bin = int(genome[chro] / ResHiC)
    for i in range(max_bin+1):
        data[i] = []
    for i in inter:
        data[i['bin1']].append((i['bin1'],i['bin2'],i['IF']))

    dtype = np.dtype({'names':['bin1','bin2','IF'],
                      'formats':[np.int, np.int, np.float]})
    for k,v in data.items():
        v = np.array(v,dtype = dtype)
        data[k] = v

    return data

def Get_inter(Index_npz,start,end):
    """
    """
#    Index_npz = read_npz(lib)
    inter = Index_npz[start]
    for i in range(end-start):
        inter = np.hstack((inter,Index_npz[start+i+1]))
    
    return inter

def Get_loops(LoopFil):
    """
    cell : String ,The name of cells  e.g.'fESC' , 'ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loops = []
    Loop = np.loadtxt(LoopFil, usecols = (0,1,2) , dtype = loop_type, skiprows = 0)
    for i in Loop:
        if i['end'] - i['start'] >= 300000:
            loops.append(i)
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops



LoopFil = sys.argv[1]
R = sys.argv[2]
Len = sys.argv[3]
point = sys.argv[4]
ave = sys.argv[5]
HiCFolder = '/public/home/xxli/data/BDF1/HiC/runHiC/workspace/Raw-GRCm38_68_chr'
mm10 = open('/public/home/xxli/data/ref/haplotype/mm10.txt' , 'r')
cells = ['CCS' , 'NT5' , 'NT6' , 'fESC']
cell_site = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
enzyme = 'MboI'
ResHiC = int(R)
res = str(ResHiC//1000) + 'K'
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19', 'X']
loop_type = np.dtype({'names':['chr','start','end'],
                     'formats':['S2', np.int, np.int]})
uniontype = np.dtype({'names':['chr','start','end','CCS','NT5','NT6','fESC'],
                          'formats':['S2' , np.int , np.int , np.float , np.float , np.float , np.float]})
#---------------------------------union_loop--------------------------------
union = Get_loops(LoopFil)
#--------------------------------Get_loop_strength--------------------------


genome ={}
for i in mm10:
    i = i.strip().split()
    genome[i[0]] = int(i[1])


Out = open(point,'w')
Out1 = open(ave,'w')
Out.writelines('\t'.join(['Chr' , 'Start' , 'End' , 'CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')
Out1.writelines('\t'.join(['Chr' , 'Start' , 'End' , 'CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')

Loop_strength = {}
Loop_strength_1 = {}
for i in union:
    Loop_strength[(i[0] , i[1] , i[2])] = [0 , 0 , 0 ,0]
    Loop_strength_1[(i[0] , i[1] , i[2])] = [0 , 0 , 0 ,0]
for c in cells:
    Pre = '-'.join([c, enzyme, 'allReps-filtered', res])
    HiCFil = Pre + '-sparse.npz'
    HiCSource = os.path.join(HiCFolder, HiCFil)
    lib = np.load(HiCSource)
    for g in chrom:
        Index_npz = read_npz(g,lib)
        for l in union[union['chr'] == g]:
            chro = l['chr']
            start = int(l['start'] / ResHiC)
            end = int(l['end'] / ResHiC)
            inter = Get_inter(Index_npz,start-10,start+10)
            (point,ave) = Get_Loop_Strength(inter,start,end,int(Len))
            Loop_strength[(l[0] , l[1] , l[2])][cell_site[c]] = point
            Loop_strength_1[(l[0] , l[1] , l[2])][cell_site[c]] = ave
            
for i in Loop_strength:
    Out.writelines('\t'.join([i[0] , str(i[1]) , str(i[2]) , str(Loop_strength[i][0]) , str(Loop_strength[i][1]) , 
                            str(Loop_strength[i][2]) , str(Loop_strength[i][3])]) + '\n')
for i in Loop_strength_1:
    Out1.writelines('\t'.join([i[0] , str(i[1]) , str(i[2]) , str(Loop_strength_1[i][0]) , str(Loop_strength_1[i][1]) , 
                            str(Loop_strength_1[i][2]) , str(Loop_strength_1[i][3])]) + '\n')
Out.close()
Out1.close()
