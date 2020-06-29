# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 16:39:31 2019

@author: han-luo
"""

from __future__ import division
import math
import numpy as np
import csv , copy
import xlrd
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import seaborn as sns
import copy
import scipy
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Our Own Color Map
#my_cmap = plt.get_cmap('bwr')
#my_cmap.set_bad('#2672a1')

#my_cmap = LinearSegmentedColormap.from_list('interaction',
#                                            ['mediumblue' , 'white' , 'red'])
#my_cmap.set_bad('#2672a1')
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['skyblue' , 'k' , 'yellow'])
my_cmap.set_bad('#2672a1')


pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['S8' , np.float]})
data_type = np.dtype({'names':['CCs' , 'NT5' , 'NT6' , 'F40'] , 
                    'formats':[np.float , np.float , np.float , np.float]})




def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix



def K_means_cluster(matrix_0,matrix,n):
    kmeans = KMeans(n_clusters=n, random_state=0).fit(matrix)
    a = kmeans.labels_
    a = a.reshape(len(a),1)

    matrix_1 = np.array(np.hstack((matrix_0,a)))
    matrix_1 = np.array(sorted(matrix_1, key=lambda x:x[-1]))
    return matrix_1

    
def plot_heatmap(matrix , vmin , vmax):
    left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
    im = ax.imshow(matrix,vmax=vmax,vmin=vmin,cmap=my_cmap,aspect = 'auto')
    x = ['CCS','NT5','NT6','fESC']
    ax.set_xticks(np.arange(len(x)))
    ax.set_xticklabels(x,fontsize = 10)
    ax.set_ylabel('c1:' + str(len(c1)) + ' , c2:' + str(len(c2)) + ' , c3:' + str(len(c3)) + ' , c4:' + str(len(c4)))
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('Compartment numbers:' + str(len(matrix)),fontsize = 30)
    ##Colorbar
    ax = fig.add_axes([left + width + 0.035 , bottom , 0.035 , 0.1])
    cbar = fig.colorbar(im,cax = ax, orientation='vertical')
    cbar.set_ticks([vmin , vmax])
    
    return fig

def rset_classify(matrix , order):
    matrix_new = matrix[matrix[:,-1] == order[0]]
    for i in order[1:]:
        m = matrix[matrix[:,-1] == i]
        matrix_new = np.array(np.vstack((matrix_new ,m)))

    return matrix_new

def Load_gtf(gtfil):
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['S64' , 'S64' , 'S8' , 'S4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf

def Sort_array(matrix , col1 , col2):
    '''
    a: list of needed to be sorted
    '''
    index = np.lexsort([matrix[:,col1] , matrix[: , col2]])
    matrix = matrix[index]
    return matrix


def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['chr' , 'start' , 'end' , 'CCS_pc1' , 'NT5_pc1' , 'NT6_pc1' , 'fESC_pc1']) + '\n')
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
   
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()

                    
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = pc_type)
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = pc_type)
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = pc_type)

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000


A_B = []
B_A = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NT5 = NT5[NT5['chr'] == g]
    tmp_NT6 = NT6[NT6['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    for i in range(len(tmp_CCS)):
        if (tmp_CCS[i]['pc'] > 0) and (tmp_NT5[i]['pc'] < 0) and (tmp_NT6[i]['pc'] < 0) and (tmp_fESC[i]['pc'] < 0):
            A_B.append([tmp_CCS[i]['pc'] , tmp_NT5[i]['pc'] , tmp_NT6[i]['pc'] , tmp_fESC[i]['pc'] , int(g) , i])
        elif (tmp_CCS[i]['pc'] < 0) and (tmp_NT5[i]['pc'] > 0) and (tmp_NT6[i]['pc'] > 0) and (tmp_fESC[i]['pc'] > 0):
            B_A.append([tmp_CCS[i]['pc'] , tmp_NT5[i]['pc'] , tmp_NT6[i]['pc'] , tmp_fESC[i]['pc'] , int(g) , i])
            
matrix = []
matrix.extend(A_B)
matrix.extend(B_A)
matrix_0 = np.array(matrix)
index_0 = matrix_0[:,4:]


matrix_1 = Z_score(matrix_0[:,:4])
index_1 = np.arange(len(matrix_1)).reshape(len(matrix_1) , 1)
matrix = np.hstack((matrix_1 , index_1))

matrix = K_means_cluster(matrix , matrix[:,:4] , 10)

#matrix = rset_classify(matrix , [3,1,6,7,2,0,4,5,8,9])


c1 = [] ; c2 = [] ; c3 = [] ; c4 = []
c1 = matrix[matrix[:,5] == 3] 
c3 = matrix[matrix[:,5] == 2]

c3 = Sort_array(c3 , 2 , 3)
c1 = Sort_array(c1 , 1 , 2)

for i in range(len(A_B)):
    if i not in c1[:, 4]:
        c2.append(list(matrix_1[i][:4]) + [i , 0])
c2 = np.array(c2)

for i in range(len(B_A)):
    if i + len(A_B) not in c3[:,4]:
        c4.append(list(matrix_1[i + len(A_B)][:4]) + [i + len(A_B) , 1])
        
c4 = np.array(c4)
matrix = np.vstack((c1 , c2 , c3 , c4))
fig = plot_heatmap(matrix[:,:4] , -1.5 , 1.5)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\D_Compartment_heatmap_Kmeans10_1.pdf')


pos1 = [] ; pos2 = [] ; pos3 = [] ; pos4 = []
for i in c1:
    index = int(i[4])
    chro = str(int(index_0[index][0]))
    start = int(index_0[index][1] * res)
    end = int((index_0[index][1] + 1) * res)
    pos1.append((chro , start , end , matrix_0[index][0] , matrix_0[index][1] , matrix_0[index][2] , matrix_0[index][3]))
    
for i in c2:
    index = int(i[4])
    chro = str(int(index_0[index][0]))
    start = int(index_0[index][1] * res)
    end = int((index_0[index][1] + 1) * res)
    pos2.append((chro , start , end , matrix_0[index][0] , matrix_0[index][1] , matrix_0[index][2] , matrix_0[index][3]))
    

for i in c3:
    index = int(i[4])
    chro = str(int(index_0[index][0]))
    start = int(index_0[index][1] * res)
    end = int((index_0[index][1] + 1) * res)
    pos3.append((chro , start , end , matrix_0[index][0] , matrix_0[index][1] , matrix_0[index][2] , matrix_0[index][3]))
    

for i in c4:
    index = int(i[4])
    chro = str(int(index_0[index][0]))
    start = int(index_0[index][1] * res)
    end = int((index_0[index][1] + 1) * res)
    pos4.append((chro , start , end , matrix_0[index][0] , matrix_0[index][1] , matrix_0[index][2] , matrix_0[index][3]))
    
Write2fils_nochr('H:\Workspace_New\data\HiC\Compartment\compartment_new\compartment_classify\\Compartment_cluster1_pc1.txt' , pos1)
Write2fils_nochr('H:\Workspace_New\data\HiC\Compartment\compartment_new\compartment_classify\\Compartment_cluster2_pc1.txt' , pos2)
Write2fils_nochr('H:\Workspace_New\data\HiC\Compartment\compartment_new\compartment_classify\\Compartment_cluster3_pc1.txt' , pos3)
Write2fils_nochr('H:\Workspace_New\data\HiC\Compartment\compartment_new\compartment_classify\\Compartment_cluster4_pc1.txt' , pos4)









