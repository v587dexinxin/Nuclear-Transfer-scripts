# -*- coding: utf-8 -*-
"""
Created on Sat Nov 03 00:48:53 2018

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

#----------------------------------------------Function--------------------------------------------------------------------------------

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a

def K_means_cluster(matrix_0,matrix,n):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1]))
    return matrix_1

def Load_loop_strength(LoopFil):
    union_type = np.dtype({'names':['chr','start','end','CCS','NT5','NT6','fESC'],
                       'formats':['S2' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_sites = np.loadtxt(LoopFil , skiprows = 1 , dtype = union_type)        
    sort_union_sites = sorted(union_sites , key = lambda x:(x[0],x[1]))
    union_sites = np.array(sort_union_sites , dtype = union_type)
    return union_sites

def Load_loop_strength_in_pairs(union_sites):
    al = {}
    ccs = []
    nt5 = []
    nt6 = []
    fesc = []
    loop = []
    for i in union_sites:
        loop.append((i['chr'] , i['start'] , i['end']))
        ccs.append(i['CCS'])
        nt5.append(i['NT5'])
        nt6.append(i['NT6'])
        fesc.append(i['fESC'])
    

    for i in range(len(ccs) - 1 , -1 , -1):
        if (ccs[i] > 50) or (nt5[i] > 50) or (nt6[i] > 50) or (fesc[i] > 50):
#        print i
            del ccs[i]
            del nt5[i]
            del nt6[i]
            del fesc[i]
            del loop[i]

    al = {'CCS_NT5':[ccs,nt5] , 'CCS_NT6':[ccs,nt6] , 'CCS_fESC':[ccs,fesc] , 'NT5_NT6':[nt5,nt6] , 'NT5_fESC':[nt5,fesc] , 'NT6_fESC':[nt6,fesc]}
    return ccs , nt5 , nt6 , fesc , loop , al

def union_loop(a,peak_type):
    '''
    a: np.array of peaks ,dtype = peak_type
    '''
    peak_new = []
    a = np.array(a , dtype = peak_type)
    for i in a:
        s = i['start']
        e = i['end']
        mask = (s == a[a['chr'] == i['chr']]['start']) & (e == a[a['chr'] == i['chr']]['end'])
        overlap = a[a['chr'] == i['chr']][mask]
        if overlap.size != 0:
            peak_new.append(i)
            a = list(a)
            for j in overlap:
                a.remove(j)
            a = np.array(a, dtype = peak_type)
        else:
            continue

    peak_new = Sort(peak_new , 0 ,1)
    peak_new = np.array(peak_new , dtype = peak_type)
    return peak_new

    
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
         
def scatter_loop_Plot(cell , data , title , fc = 1):
    datatype = ({'names':['l1','l2','FC'],
                 'formats':[np.float,np.float,np.float]})
    loop = []
    for i in range(len(data[0])):
        l1 = data[0][i]
        l2 = data[1][i]
        
        p = np.log2((l2+0.00001)/(l1+0.00001))
        loop.append((l1,l2,p))
    loop = np.array(loop , dtype = datatype)
    r = np.corrcoef(loop['l1'] ,loop['l2'])
    l1_bound = fc
    l2_bound = -fc
    
    l1_mask = ((loop['l1'] >= 3) | (loop['l2'] >= 3) ) & (loop['FC'] <= l2_bound)
#    l1_mask = (loop['FC'] <= l2_bound)
    l1_loops = loop[l1_mask]
    
    l2_mask = ((loop['l1'] >= 3) | (loop['l2'] >= 3) ) & (loop['FC'] >= l1_bound)
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
    ax.plot([0 , 3] , [3 , 3] , ls = '--' , c = 'gray' , lw = 1.0 )
    ax.plot([3 , 3] , [0 , 3] , ls = '--' , c = 'gray' , lw = 1.0 )
    
    ax.set_xlabel(cell.split("_")[0] + '_loop_strength' , size = 23)
    ax.set_ylabel(cell.split("_")[1] + '_loop_strength' , size = 23)
    ax.set_xlim(0,30)
    ax.set_ylim(0,30)
    ax.text(25 , 3 , cell.split("_")[0] + ': %d' % len(l1_loops) , size = 20)
    ax.text(3 , 25 , cell.split("_")[1] + ': %d' % len(l2_loops) , size = 20)
    ax.text(3 , 26 , 'r = %f' % r[0][1] , size = 20)
    ax.set_title(title , size = 25)
    return fig

def Heat_map(matrix , title):            
    left, bottom, width, height = 0.1, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    #matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    #matrix = matrix_1[:,:-1]
    vmax = np.percentile(matrix,90)
    vmin = np.percentile(matrix,10)
    im = ax.imshow(matrix , vmax = vmax , vmin = vmin , cmap=my_cmap,aspect = 'auto')
    plt.colorbar(im)
    x = ['CCS','NT5','NT6','fESC']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 15)
    plt.title(title,fontsize = 20)



#---------------------------------------------scatter-----------------------------------------------------------------------------------
LoopFil = 'D:\\workspace\\loop\\loop_strength_20K\\union_sites_ave_20K.txt'
union_sites = Load_loop_strength(LoopFil)
(ccs , nt5 , nt6 , fesc , loop , al) = Load_loop_strength_in_pairs(union_sites)
for c in al:
    fig1 = scatter_loop_Plot(c , al[c] , c.split("_")[0] + ' & ' + c.split("_")[1] , 1)
    plt.savefig('D:\\workspace\\loop\\scatter\\loop_strength_new\\loop_strength_20K_union\\' + c + '_loop_strength_ave_20K.png')
    


#-------------------------------------------matrix-------------------------------------------------------------------------

diff_loop = Select_diff_loop(loop , al)
    

matrix = np.zeros((len(diff_loop) , 4))
for i in range(len(matrix)):
    index = loop.index((diff_loop[i][0] , diff_loop[i][1] , diff_loop[i][2]))
    matrix[i][0] = ccs[index]
    matrix[i][1] = nt5[index]
    matrix[i][2] = nt6[index]
    matrix[i][3] = fesc[index]
    
##--------------------------------------------Z-score 标准化------------------------------------------------------------------------   
matrix = stats.zscore(matrix , axis = 1 , ddof = 1)
matrix[np.isnan(matrix)] = 0

matrix_1 = K_means_cluster(matrix , matrix , 8)  

#------------------------------------------Heatmap_Plot----------------------------------------------------------------------------

Heat_map(matrix_1[:,:-1] , 'diff_loop_strength loop numbers:' + str(len(matrix)))