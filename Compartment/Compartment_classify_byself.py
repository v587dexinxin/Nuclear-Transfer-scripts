# -*- coding: utf-8 -*-
"""
Created on Mon Jan 06 22:22:40 2020

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
                    
 

def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a

def Z_score(matrix):
    matrix = stats.zscore(matrix, axis=1, ddof=1)
    matrix[np.isnan(matrix)] = 0 
    return matrix

                   
def Classify4_compartment(comp1 , comp2 , comp3 , comp4):
    classify = {'AAA':[] , 'AAB':[] , 'ABA':[] , 'ABB':[] , 'BBB':[] , 'BBA':[] , 'BAA':[] , 'BAB':[]}
    for g in chrom:
        tmp_comp1 = comp1[comp1['chr'] == g]
        tmp_comp2 = comp2[comp2['chr'] == g]
        tmp_comp3 = comp3[comp3['chr'] == g]
        tmp_comp4 = comp4[comp4['chr'] == g] 
        for i in range(len(tmp_comp1)):
            score1 = tmp_comp1[i]['pc']
            score2 = tmp_comp2[i]['pc']
            score3 = tmp_comp3[i]['pc']
            score4 = tmp_comp4[i]['pc']
            if (score1 > 0) and (score2 > 0) and (score3 > 0) and (score4 > 0):
                classify['AAA'].append((score1 , score2 , score3 , score4, g , i))
            elif (score1 > 0) and (score2 > 0) and (score3 > 0) and (score4 < 0):
                classify['AAB'].append((score1 , score2 , score3 , score4 , g , i))
            elif (score1 > 0) and (score2 < 0) and (score3 < 0) and (score4 > 0):
                classify['ABA'].append((score1 , score2 , score3 , score4 , g , i))
            elif (score1 > 0) and (score2 < 0) and (score3 < 0) and (score4 < 0):
                classify['ABB'].append((score1 , score2 , score3 , score4 ,g , i))
            elif (score1 < 0) and (score2 < 0) and (score3 < 0) and (score4 < 0):
                classify['BBB'].append((score1 , score2 , score3 , score4 , g , i))
            elif (score1 < 0) and (score2 < 0) and (score3 < 0) and (score4 > 0):
                classify['BBA'].append((score1 , score2 , score3 , score4 , g , i))
            elif (score1 < 0) and (score2 > 0) and (score3 > 0) and (score4 > 0):
                classify['BAA'].append((score1 , score2 , score3 , score4 , g , i))
            elif (score1 < 0) and (score2 > 0) and (score3 > 0) and (score4 < 0):
                classify['BAB'].append((score1 , score2 , score3 , score4 , g , i))
    return classify


def get_4cell_matrix(matrix_raw):
    matrix = np.zeros((len(matrix_raw) , 4))
    pos = []
    n = 0
    pos.append(('chr' , 'start' , 'end' , 'CCS_pc1' , 'NT5_pc1' , 'NT6_pc1' , 'fESC_pc1'))
    for i in matrix_raw:
        g = i[-2]
        start = i[-1] * 200000
        end = (i[-1] + 1) * 200000
        matrix[n , 0] = i[0]
        matrix[n , 1] = i[1]
        matrix[n , 2] = i[2]
        matrix[n , 3] = i[3]
        pos.append((g , start , end , i[0] , i[1] , i[2] , i[3]))
        n += 1
    return matrix , pos
    

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
    ax.set_ylabel('c1:' + str(len(Matrix1)) + ' , c2:' + str(len(Matrix2)) + ' , c3:' + str(len(Matrix3)) + ' , c4:' + str(len(Matrix4)) + ' , c5:' + str(len(Matrix5)) + ' , c6:' + str(len(Matrix6)) + ' , c7:' + str(len(Matrix7)) + ' , c8:' + str(len(Matrix8)))
    ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
    plt.title('Compartment numbers:' + str(len(matrix)),fontsize = 30)
    ##Colorbar
    ax = fig.add_axes([left + width + 0.035 , bottom , 0.035 , 0.1])
    cbar = fig.colorbar(im,cax = ax, orientation='vertical')
    cbar.set_ticks([vmin , vmax])
    
    return fig



def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()

        
    
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = pc_type)
NT5 = NT5[NT5['chr'] != 'X']
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = pc_type)
NT6 = NT6[NT6['chr'] != 'X']
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']



chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']

##mean_NT5_NT6



classify = Classify4_compartment(CCS , NT5 , NT6 , fESC)


matrix1 = [] ; matrix2 = [] ; matrix3 = [] ; matrix4 = [] ; matrix5 = [] ; matrix6 = [] ; matrix7 = [] ; matrix8 = []

Sort(classify['ABB'] , 0 , 0) ; Sort(classify['BAA'] , 0 , 0)
for i in classify['ABB']:
    if (i[1] + i[2]) / 2 > 0.5 * i[3]:
        matrix1.append(i)
    elif (i[1] + i[2]) / 2 < 2 * i[3]:
        matrix2.append(i)
    else:
        matrix3.append(i)
        
matrix4 = Sort(classify['AAB'] , 0 , 0)        

for i in classify['BAA']:
    if (i[1] + i[2]) / 2 < 0.5 * i[3]:
        matrix5.append(i)
    elif (i[1] + i[2]) / 2 > 2 * i[3]:
        matrix6.append(i)
    else:
        matrix7.append(i)
        
matrix8 = Sort(classify['BBA'] , 0 , 0)
        
Matrix1,pos1 = get_4cell_matrix(matrix3) ; Matrix2,pos2 = get_4cell_matrix(matrix1)
Matrix3,pos3 = get_4cell_matrix(matrix2) ; Matrix4,pos4 = get_4cell_matrix(matrix4) 
Matrix5,pos5 = get_4cell_matrix(matrix7) ; Matrix6,pos6 = get_4cell_matrix(matrix5)
Matrix7,pos7 = get_4cell_matrix(matrix6) ; Matrix8,pos8 = get_4cell_matrix(matrix8)
 
#pos3 = [pos3[0]] + pos3[6:] + pos3[1:6]

Matrix1 = Z_score(Matrix1) ; Matrix2 = Z_score(Matrix2) ; Matrix3 = Z_score(Matrix3) 
Matrix4 = Z_score(Matrix4) ; Matrix5 = Z_score(Matrix5) ; Matrix6 = Z_score(Matrix6) 
Matrix7 = Z_score(Matrix7) ; Matrix8 = Z_score(Matrix8)

tmp1 = [(i[0] , i[1] , i[2] , i[3]) for i in Matrix1]
tmp2 = [(i[0] , i[1] , i[2] , i[3]) for i in Matrix2]
tmp3 = [(i[0] , i[1] , i[2] , i[3]) for i in Matrix3]
tmp4 = [(i[0] , i[1] , i[2] , i[3]) for i in Matrix4]
tmp5 = [(i[0] , i[1] , i[2] , i[3]) for i in Matrix5]
tmp6 = [(i[0] , i[1] , i[2] , i[3]) for i in Matrix6]
tmp7 = [(i[0] , i[1] , i[2] , i[3]) for i in Matrix7]
tmp8 = [(i[0] , i[1] , i[2] , i[3]) for i in Matrix8]

 
Matrix1 = Sort(tmp1 , 0 , 3) ; Matrix2 = Sort(tmp2 , 0 , 3) ; Matrix3 = Sort(tmp3 , 0 , 3)
Matrix4 = Sort(tmp4 , 0 , 3) ; Matrix5 = Sort(tmp5 , 0 , 3) ; Matrix6 = Sort(tmp6 , 0 , 3)
Matrix7 = Sort(tmp7 , 0 , 3) ; Matrix8 = Sort(tmp8 , 0 , 3)
#Matrix3 = np.vstack((Matrix3[5:] , Matrix3[:5])) 

matrix = np.vstack((Matrix1 , Matrix2 , Matrix3 , Matrix4 , Matrix5 , Matrix6 , Matrix7 , Matrix8))

#matrix = Z_score(matrix)
fig =plot_heatmap(matrix , -1, 1)

run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\D_Compartment_heatmap_byself.pdf')
          
Write2fils_nochr('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster1_pc1.txt' , pos1)
Write2fils_nochr('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster2_pc1.txt' , pos2)
Write2fils_nochr('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster3_pc1.txt' , pos3)
Write2fils_nochr('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster4_pc1.txt' , pos4)
Write2fils_nochr('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster5_pc1.txt' , pos5)
Write2fils_nochr('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster6_pc1.txt' , pos6)
Write2fils_nochr('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster7_pc1.txt' , pos7)
Write2fils_nochr('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\compartment_classify_byself_new\\Compartment_cluster8_pc1.txt' , pos8)


#-------------------------plot_compartment_states-------------------------------
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['mediumblue' , 'orange'])
my_cmap.set_bad('#2672a1')

matrix = np.array([[1 , -1 , -1] , [1 , -1 , -1] , [1 , -1 , -1] , [1 , 1 , -1] , [-1 , 1 , 1] , [-1 , 1 , 1] , [-1 , 1 , 1] , [-1 , -1 , 1]])
left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (6, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))
im = ax.imshow(matrix,vmax=1,vmin=-1,cmap=my_cmap,aspect = 'auto')
x = ['CCS','NT','fESC']
ax.set_xticks(np.arange(len(x)))
ax.set_xticklabels(x,fontsize = 10)
ax.set_yticks(np.arange(len(matrix)))
ax.set_yticklabels(['ABrB' , 'ABpB' , 'ABoB' , 'AAB' , 'BArA' , 'BApA' , 'BAoA' , 'BBA'],fontsize = 10)
ax.tick_params(axis = 'y', labelsize = 18, pad = 7)
for i in range(len(matrix)):
    ax.plot([0 - 0.5 , 3 - 0.5] , [i - 0.5 , i - 0.5] , c = 'black')
for j in range(len(x)):
    ax.plot([j - 0.5 , j - 0.5] , [0 - 0.5 , 8 - 0.5] , c = 'black')
        
plt.title('Compartment states',fontsize = 30)

    
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\D_Compartment_states.pdf')
      
 
#-------------------------plot_compartment_states_number-------------------------------               
c1 = len(Matrix1) / 546 ; c2 = len(Matrix2) / 546 ; c3 = len(Matrix3) / 546 ; c4 = len(Matrix4) / 546
c5 = len(Matrix5) / 1889 ; c6 = len(Matrix6) / 1889 ; c7 = len(Matrix7) / 1889 ; c8 = len(Matrix8) / 1889

left, bottom, width, height = 0.2, 0.1, 0.60, 0.80
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (6, 12))
ax = fig.add_axes(size_axes)
#matrix_1 = np.hstack((matrix_b[:,:-1] , gene_loop_matrix_1))

ax.bar(1 , c1 + c2 + c3 + c4 , color='green')
ax.bar(1 , c1 + c2 + c3 , color='red')
ax.bar(1 , c1 + c2 , color='orange')
ax.bar(1 , c1 , color='mediumblue')

ax.bar(2 , c5 + c6 + c7 + c8 , color='green')
ax.bar(2 , c5 + c6 + c7 , color='red')
ax.bar(2 , c5 + c6 , color='orange')
ax.bar(2 , c5 , color='mediumblue')



x = ['A-B','B-A']
ax.set_xticks([1 , 2])
ax.set_xticklabels(x,fontsize = 10)
ax.set_ylabel('Percent(%)',fontsize = 10)
ax.set_ylim((0.5 , 2.5))
ax.set_ylim((0 , 1))

        
plt.title('Compartment states numbers',fontsize = 10)

    
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig2_figs\\D_Compartment_states_numbers.pdf') 
