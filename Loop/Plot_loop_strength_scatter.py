# -*- coding: utf-8 -*-
"""
Created on Fri Dec 07 16:18:08 2018

@author: xxli
"""

from __future__ import division
import numpy as np
import sys, cPickle
import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math
from sklearn.cluster import KMeans
from scipy import stats

## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'


# Our Own Color Map
my_cmap = plt.get_cmap('bwr')
my_cmap.set_bad('#D3D3D3')


def Load_loop_strength(LoopFil):
    union_type = np.dtype({'names':['chr','start','end','CCS','NT5','NT6','fESC'],
                       'formats':['S2' , np.int , np.int , np.float , np.float , np.float , np.float]})
    
    union_sites = np.loadtxt(LoopFil , skiprows = 1 , dtype = union_type)        
    new = []
    for i in union_sites:
        if i['end'] - i['start'] >= 300000:
            new.append(i)
        else:
            print i
            continue
    
    sort_union_sites = sorted(new , key = lambda x:(x[0],x[1]))
    union_sites = np.array(sort_union_sites , dtype = union_type)
    return union_sites


def Select_diff_loop(loop , al):
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['S8' , np.int , np.int]})
    vs_loop = {}    
    for c in al:
        vs_loop[c] = []
        for i in range(len(al[c][0])):
            l1 = al[c][0][i]
            l2 = al[c][1][i]
            vs_loop[c].append((loop[i][0] , float(loop[i][1]) , float(loop[i][2]) , l1 , l2 ))
        
    diff_loop = []
    for c in vs_loop:
        for i in vs_loop[c]:
            if (i[-1] >= 3) | (i[-2] >= 3) and ((i[-1] >= 2 * i[-2]) or (i[-1] <= 0.5 * i[-2])):
                diff_loop.append((i[0] , float(i[1]) , float(i[2])))
            
    diff_loop = union_loop(diff_loop , loop_type)   
    return diff_loop
         
def scatter_loop_Plot(cell , data , fc = 1):
    datatype = ({'names':['l1','l2','log2FC'],
                 'formats':[np.float,np.float,np.float]})
    loop = []
    for i in range(len(data[0])):
        l1 = data[0][i]
        l2 = data[1][i]
        
        p = (l2+0.00001)/(l1+0.00001)
        loop.append((l1,l2,p))
    loop = np.array(loop , dtype = datatype)
    r = np.corrcoef(loop['l1'] ,loop['l2'])
    l1_bound = fc
    l2_bound = 1/fc
    
    l1_mask = ((loop['l1'] >= 0.58496) | (loop['l2'] >= 0.58496) ) & (loop['log2FC'] <= l2_bound)
#    l1_mask = (loop['FC'] <= l2_bound)
    l1_loops = loop[l1_mask]
    
    l2_mask = ((loop['l1'] >= 0.58496) | (loop['l2'] >= 0.58496) ) & (loop['log2FC'] >= l1_bound)
#    l2_mask = (loop['FC'] >= l1_bound)
    l2_loops = loop[l2_mask]
    
    Non_loop = loop[~(l1_mask | l2_mask)]
    
    left, bottom, width, height = 0.1, 0.1, 0.80, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_axes(size_axes)
    ax.scatter(l1_loops['l1'] , l1_loops['l2'] , alpha = 0.8 , c = 'red')
    ax.scatter(l2_loops['l1'] , l2_loops['l2'] , alpha = 0.8 , c = 'blue')
    ax.scatter(Non_loop['l1'] , Non_loop['l2'] , alpha = 0.5 , c = 'gray')
    ax.plot([0 , 50] , [0 , 50] , ls = '--' , c = 'gray' , lw = 1.0 )
    ax.plot([0 , 50] , [0 , 100] , ls = '--' , c = 'gray' , lw = 1.0 )
    ax.plot([0 , 50] , [0 , 25] , ls = '--' , c = 'gray' , lw = 1.0 )
    ax.plot([0 , 1.5] , [1.5 , 1.5] , ls = '--' , c = 'gray' , lw = 1.0 )
    ax.plot([1.5 , 1.5] , [0 , 1.5] , ls = '--' , c = 'gray' , lw = 1.0 )
    
    ax.set_xlabel(cell.split("_")[0] + '_loop_strength' , size = 23)
    ax.set_ylabel(cell.split("_")[1] + '_loop_strength' , size = 23)
#    ax.set_xticks([0 , 5 , 10 , 15 , 20 , 25 , 30])
#    ax.set_xticks([5 , 10 , 15 , 20 , 25 , 30])
    ax.set_xlim(0,15)
    ax.set_ylim(0,15)
    ax.text(11 , 2 , cell.split("_")[0] + ': %d' % len(l1_loops) , size = 20)
    ax.text(2 , 13 , cell.split("_")[1] + ': %d' % len(l2_loops) , size = 20)
    ax.text(2 , 14 , 'r = %f' % r[0][1] , size = 20)
    ax.set_title(c.split('_')[0] + '_vs_' + c.split("_")[1] + ' Numbers:' + str(len(data[0]))  , size = 25)
    return fig


cell = ['CCS_NT5' , 'CCS_NT6' , 'CCS_fESC' , 'NT5_NT6' , 'NT5_fESC' , 'NT6_fESC']
al = {}
for c in cell:
    c1 = c.split('_')[0] ; c2 = c.split('_')[1]
    r1 = [] ; r2 = []
    loop_strength = Load_loop_strength('D:/Workspace_New/data/HiC/Loop/union_inpairs_loop_strength/' + c + '_ave.txt')
    for i in loop_strength:
        if (((i[c1] >= 1.5) and (i[c1] <= 10)) or ((i[c2] >= 1.5) and (i[c2] <= 10))) : 
            r1.append(i[c1])
            r2.append(i[c2])
    al[c] = [r1 , r2]
    
    
pp = PdfPages('D:/Workspace_New/Plot/Figure1/Correlation/diff_cells_scatter/loops/loop_strength_ave.pdf')    
for c in cell:
    fig = scatter_loop_Plot(c , al[c] , fc = 2)
    pp.savefig(fig)
pp.close()
    