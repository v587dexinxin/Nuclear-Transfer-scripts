# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 20:11:29 2021

@author: xxli
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
from matplotlib_venn import venn2
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
                    'formats':['U8' , np.float]})
data_type_1 = np.dtype({'names':['chr' , 'pos' , 'CCs' , 'NTs' , 'fESC'] , 
                        'formats':['U8' , np.int , np.float , np.float , np.float]})
data_type_2 = np.dtype({'names':['chr' , 'pos' , 'MEF' , 'IPSC' , 'E14'] , 
                        'formats':['U8' , np.int , np.float , np.float , np.float]})
pc_type_1 = np.dtype({'names':['chr' , 'pos' , 'pc'] , 
                        'formats':['U8' , np.int , np.float]})



def get_PC(PC_fil):
    '''
    '''
    chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
    PC_new = []
    PC = np.loadtxt(PC_fil , dtype = pc_type)
    PC = PC[PC['chr'] != 'X']
    for g in chroms:
        tmp = PC[PC['chr'] == g]
        for i in range(len(tmp)):
            PC_new.append((g , i * 200000 , tmp[i]['pc']))
    PC_new = np.array(PC_new , dtype = pc_type_1)
    
    return PC_new

    
def get_AB_percent(data):
    data_1 = {'CCS' : [0 , 0] , 'NTs' : [0 , 0] , 'fESC' : [0 , 0]}
    for i in data:
        if i[2] > 0:
            data_1['CCS'][0] += 1
        elif i[2] < 0:
            data_1['CCS'][1] += 1
        
        if i[3] > 0:
            data_1['NTs'][0] += 1
        elif i[3] < 0:
            data_1['NTs'][1] += 1
            
        if i[4] > 0:
            data_1['fESC'][0] += 1
        elif i[4] < 0:
            data_1['fESC'][1] += 1
            
    A_percent = {}
    for c in data_1:
        A_percent[c] = data_1[c][0] / (data_1[c][0] + data_1[c][1])
    return A_percent
    
    
    
def get_IPSC_AB_percent(data):
    data_1 = {'MEF' : [0 , 0] , 'IPSC' : [0 , 0] , 'E14' : [0 , 0]}
    for i in data:
        if i[2] > 0:
            data_1['MEF'][0] += 1
        elif i[2] < 0:
            data_1['MEF'][1] += 1
        
        if i[3] > 0:
            data_1['IPSC'][0] += 1
        elif i[3] < 0:
            data_1['IPSC'][1] += 1
            
        if i[4] > 0:
            data_1['E14'][0] += 1
        elif i[4] < 0:
            data_1['E14'][1] += 1
            

    A_percent = {}
    for c in data_1:
        A_percent[c] = data_1[c][0] / (data_1[c][0] + data_1[c][1])
    return A_percent
    
def Bar_plot_NT(x , yA , yB , label):
    left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    if label == 'A_B': 
        ax.bar(x , yB , color = 'darkblue')
        ax.bar(x , yA , color = 'darkorange')
    else:
        yA = [1 - i for i in yA]   
        ax.bar(x , yB , color = 'darkorange')
        ax.bar(x , yA , color = 'darkblue')
        
    xticks = x
    labels = ['CCS' , 'NTs' , 'fESC']
    for a, b in zip(x, yA):
        plt.text(a, b + 0.05, '%.4f' % b, ha='center', va='bottom', fontsize=20)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels,fontsize = 10)
    ax.set_xlabel('Cell type' , fontsize = 20 )
    ax.set_xlim((-0.5, 2.5))
    ax.set_ylim((0 , 1.01))
    return fig
    
def Bar_plot_IPS(x , yA , yB , label):
    left, bottom, width, height = 0.2, 0.2, 0.6, 0.6
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    if label == 'A_B':
        ax.bar(x , yB , color = 'darkblue')
        ax.bar(x , yA , color = 'darkorange')
    else:
        yA = [1 - i for i in yA]   
        ax.bar(x , yB , color = 'darkorange')
        ax.bar(x , yA , color = 'darkblue')        
    xticks = x
    labels = ['MEF' , 'IPS_P3' , 'E14']
    for a, b in zip(x, yA):
        plt.text(a, b + 0.05, '%.4f' % b, ha='center', va='bottom', fontsize=20)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels,fontsize = 10)
    ax.set_xlabel('Cell type' , fontsize = 20 )
    ax.set_xlim((-0.5, 2.5))
    ax.set_ylim((0 , 1.01))
    return fig        
   
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()


def Sort(a , s_1):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1]))
    return a


                    
CCS = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\CCS_Traditonal_PC_200K_Compartment_200K.txt')
NT5 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NT5_Traditonal_PC_200K_Compartment_200K.txt')
NT6 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NT6_Traditonal_PC_200K_Compartment_200K.txt')
F35 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\F35_Traditonal_PC_200K_Compartment_200K.txt')
F40 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\F40_Traditonal_PC_200K_Compartment_200K.txt')
NTs = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NTs_Traditonal_PC_200K_Compartment_200K.txt')
fESC = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\fESC_Traditonal_PC_200K_Compartment_200K.txt')




MEF = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt')
IPS_P3 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt')
IPS_P20 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P20_compartment_200K.txt')
E14 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt')



chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000


fESC_E14_common = []

for g in chroms:
    tmp_fESC = fESC[fESC['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    
    for i in range(len(tmp_fESC)):
        if tmp_fESC[i]['pos'] != tmp_E14[i]['pos']:
            print (g , i)
        if (tmp_fESC[i]['pc'] * tmp_E14[i]['pc'] > 0 ) or (tmp_fESC[i]['pc'] == tmp_E14[i]['pc'] == 0):
            fESC_E14_common.append((g , i * 200000 , (tmp_fESC[i]['pc'] + tmp_E14[i]['pc']) / 2))
fESC_E14_common = np.array(fESC_E14_common , dtype = pc_type_1)

##fESC_E14_overlap_Venn2
n1 = len(fESC) - len(fESC_E14_common)
n2 = len(E14) - len(fESC_E14_common)
n3 = len(fESC_E14_common)
fig = plt.figure(figsize = (12, 12))
venn2(subsets = (n1, n2, n3), set_labels = ('fESC', 'E14'))  

IPS_P3_P20_common = []

for g in chroms:
    tmp_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_P20 = IPS_P20[IPS_P20['chr'] == g]
    
    for i in range(len(tmp_P3)):
        if tmp_P3[i]['pos'] != tmp_P20[i]['pos']:
            print (g , i)
        if (tmp_P3[i]['pc'] * tmp_P20[i]['pc'] > 0 ) or (tmp_P3[i]['pc'] == tmp_P20[i]['pc'] == 0):
            IPS_P3_P20_common.append((g , i * 200000 , (tmp_P3[i]['pc'] + tmp_P20[i]['pc']) / 2))
IPS_P3_P20_common = np.array(IPS_P3_P20_common , dtype = pc_type_1)


##IPS_P3_IPS_P20_overlap_Venn2
n1 = len(IPS_P3) - len(IPS_P3_P20_common)
n2 = len(IPS_P20) - len(IPS_P3_P20_common)
n3 = len(IPS_P3_P20_common)
fig = plt.figure(figsize = (12, 12))
venn2(subsets = (n1, n2, n3), set_labels = ('IPSC_P3', 'IPSC_P20'))  






##NT_process
A_B = []
B_A = []


for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NTs = NTs[NTs['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    tmp_common = fESC_E14_common[fESC_E14_common['chr'] == g]
    for i in range(len(tmp_CCS)):
        if (tmp_CCS[i]['pc'] > 0) and (tmp_fESC[i]['pc'] < 0):
            if i * 200000 in tmp_common['pos']:
                A_B.append((g , i * 200000 , tmp_CCS[i]['pc'] , tmp_NTs[i]['pc'] , tmp_fESC[i]['pc']))
        elif (tmp_CCS[i]['pc'] < 0) and (tmp_fESC[i]['pc'] > 0):
            if i * 200000 in tmp_common['pos']:
                B_A.append((g , i * 200000 , tmp_CCS[i]['pc'] , tmp_NTs[i]['pc'] , tmp_fESC[i]['pc']))
A_B = np.array(A_B , dtype = data_type_1)
B_A = np.array(B_A , dtype = data_type_1)



A_B_1 = get_AB_percent(A_B)
B_A_1 = get_AB_percent(B_A)


x = range(3)
yA1 = [A_B_1[c] for c in ['CCS' , 'NTs' , 'fESC']]   
yB1 = [1 , 1 , 1]   
fig = Bar_plot_NT(x , yA1 , yB1 , 'A_B')
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_percentage\\A_NT_compartment_A_B_percentage_merged.pdf')


x = range(3)
yA2 = [B_A_1[c] for c in ['CCS' , 'NTs' , 'fESC']]   
yB2 = [1 , 1 , 1 ]   
fig = Bar_plot_NT(x , yA2 , yB2 , 'B_A')
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_percentage\\A_NT_compartment_B_A_percentage_merged.pdf')



##IPSC_P3
IPSC_A_B = []
IPSC_B_A = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    tmp_common = fESC_E14_common[fESC_E14_common['chr'] == g]
    for i in range(len(tmp_MEF)):
        if (tmp_MEF[i]['pc'] > 0) and (tmp_E14[i]['pc'] < 0):
            if i * 200000 in tmp_common['pos']:
                IPSC_A_B.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))
        elif (tmp_MEF[i]['pc'] < 0) and (tmp_E14[i]['pc'] > 0):
            if i * 200000 in tmp_common['pos']:
                IPSC_B_A.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))
IPSC_A_B = np.array(IPSC_A_B , dtype = data_type_2)
IPSC_B_A = np.array(IPSC_B_A , dtype = data_type_2)
            


IPSC_A_B_1 = get_IPSC_AB_percent(IPSC_A_B)
IPSC_B_A_1 = get_IPSC_AB_percent(IPSC_B_A)

x = range(3)
yA1 = [IPSC_A_B_1[c] for c in ['MEF' , 'IPSC' , 'E14']]   
yB1 = [1 , 1 , 1]   
fig = Bar_plot_IPS(x , yA1 , yB1 , 'A_B')
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_percentage\\IPSC_P3\\A_IPSC_compartment_A_B_percentage.pdf')

x = range(3)
yA2 = [IPSC_B_A_1[c] for c in ['MEF' , 'IPSC' , 'E14']]   
yB2 = [1 , 1 , 1]   
fig = Bar_plot_IPS(x , yA2 , yB2 , 'B_A')
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_percentage\\IPSC_P3\\A_IPSC_compartment_B_A_percentage.pdf')


##IPSC_P20
IPSC_A_B = []
IPSC_B_A = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P20 = IPS_P20[IPS_P20['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    tmp_common = fESC_E14_common[fESC_E14_common['chr'] == g]
    for i in range(len(tmp_MEF)):
        if (tmp_MEF[i]['pc'] > 0) and (tmp_E14[i]['pc'] < 0):
            if i * 200000 in tmp_common['pos']:
                IPSC_A_B.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P20[i]['pc'] , tmp_E14[i]['pc']))
        elif (tmp_MEF[i]['pc'] < 0) and (tmp_E14[i]['pc'] > 0):
            if i * 200000 in tmp_common['pos']:
                IPSC_B_A.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P20[i]['pc'] , tmp_E14[i]['pc']))
IPSC_A_B = np.array(IPSC_A_B , dtype = data_type_2)
IPSC_B_A = np.array(IPSC_B_A , dtype = data_type_2)
            


IPSC_A_B_1 = get_IPSC_AB_percent(IPSC_A_B)
IPSC_B_A_1 = get_IPSC_AB_percent(IPSC_B_A)

x = range(3)
yA1 = [IPSC_A_B_1[c] for c in ['MEF' , 'IPSC' , 'E14']]   
yB1 = [1 , 1 , 1]   
fig = Bar_plot_IPS(x , yA1 , yB1 , 'A_B')
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_percentage\\IPSC_P20\\A_IPSC_compartment_A_B_percentage.pdf')

x = range(3)
yA2 = [IPSC_B_A_1[c] for c in ['MEF' , 'IPSC' , 'E14']]   
yB2 = [1 , 1 , 1]   
fig = Bar_plot_IPS(x , yA2 , yB2 , 'B_A')
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_percentage\\IPSC_P20\\A_IPSC_compartment_B_A_percentage.pdf')


##IPSC_P3_P20_common
IPSC_A_B = []
IPSC_B_A = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3_P20_common = IPS_P3_P20_common[IPS_P3_P20_common['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    tmp_common = fESC_E14_common[fESC_E14_common['chr'] == g]
    for i in range(len(tmp_MEF)):
        if (tmp_MEF[i]['pc'] > 0) and (tmp_E14[i]['pc'] < 0):
            if (i * 200000 in tmp_common['pos']) and (i * 200000 in tmp_IPS_P3_P20_common['pos']):
                IPSC_A_B.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P3_P20_common[tmp_IPS_P3_P20_common['pos'] == i * 200000][0]['pc'] , tmp_E14[i]['pc']))
        elif (tmp_MEF[i]['pc'] < 0) and (tmp_E14[i]['pc'] > 0):
            if i * 200000 in tmp_common['pos'] and (i * 200000 in tmp_IPS_P3_P20_common['pos']):
                IPSC_B_A.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P3_P20_common[tmp_IPS_P3_P20_common['pos'] == i * 200000][0]['pc'] , tmp_E14[i]['pc']))
IPSC_A_B = np.array(IPSC_A_B , dtype = data_type_2)
IPSC_B_A = np.array(IPSC_B_A , dtype = data_type_2)
            


IPSC_A_B_1 = get_IPSC_AB_percent(IPSC_A_B)
IPSC_B_A_1 = get_IPSC_AB_percent(IPSC_B_A)

x = range(3)
yA1 = [IPSC_A_B_1[c] for c in ['MEF' , 'IPSC' , 'E14']]   
yB1 = [1 , 1 , 1]   
fig = Bar_plot_IPS(x , yA1 , yB1 , 'A_B')
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_percentage\\IPSC_P3_P20_common\\A_IPSC_compartment_A_B_percentage.pdf')

x = range(3)
yA2 = [IPSC_B_A_1[c] for c in ['MEF' , 'IPSC' , 'E14']]   
yB2 = [1 , 1 , 1]   
fig = Bar_plot_IPS(x , yA2 , yB2 , 'B_A')
run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_percentage\\IPSC_P3_P20_common\\A_IPSC_compartment_B_A_percentage.pdf')





#-----------------------------------Venn2_Plot_P3---------------------------------------

data_type = np.dtype({'names':['chr' , 'pos' ,'CCS' , 'NTs' , 'fESC' , 'MEF' , 'IPSC_P3' , 'E14'] , 
                      'formats':['U8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})


union_pc = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NTs = NTs[NTs['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    for i in range(len(tmp_CCS)):
        union_pc.append((g , i * 200000 , tmp_CCS[i]['pc'] , tmp_NTs[i]['pc'] , tmp_fESC[i]['pc'] , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))

union_pc = np.array(union_pc , dtype = data_type)
            




def Common_pc_Venn2(union_pc , lab1 , lab2 , process):
    '''
    '''
    n1 = 0 ; n2 = 0 ; n3 = 0
    
    for i in union_pc:
        if i['CCS'] != 0 and i['NTs'] != 0 and i['fESC'] != 0 and i['MEF'] != 0 and i['IPSC_P3'] != 0 and i['E14'] != 0:
            if process == 'A_B':
                if i['CCS'] > 0 and i['NTs'] < 0 and i['fESC'] < 0 and i['MEF'] > 0 and i['IPSC_P3'] < 0 and i['E14'] < 0:
                    n2 += 1
                elif i['CCS'] > 0 and i['NTs'] < 0 and i['fESC'] < 0 and i['E14'] < 0:
                    n1 += 1
                elif i['MEF'] > 0 and i['IPSC_P3'] < 0 and i['E14'] < 0 and i['fESC'] < 0:
                    n3 += 1
            else:
                if i['CCS'] < 0 and i['NTs'] > 0 and i['fESC'] > 0 and i['MEF'] < 0 and i['IPSC_P3'] > 0 and i['E14'] > 0:
                    n2 += 1
                elif i['CCS'] < 0 and i['NTs'] > 0 and i['fESC'] > 0 and i['E14'] > 0:
                    n1 += 1
                elif i['MEF'] < 0 and i['IPSC_P3'] > 0 and i['E14'] > 0 and i['fESC'] > 0:
                    n3 += 1
                
                

            
    fig = plt.figure(figsize = (12, 12))
    venn2(subsets = (n1, n3, n2), set_labels = (lab1, lab2))  
    
    return fig




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

    
fig1 = Common_pc_Venn2(union_pc , 'A_B_NTs' , 'A_B_IPSC' , 'A_B')
fig2 = Common_pc_Venn2(union_pc , 'B_A_NTs' , 'B_A_IPSC' , 'B_A')

run_Plot(fig1 , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_Venn2\\IPSC_P3\\A_B_compartment_switch_Venn2.pdf')
run_Plot(fig2 , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_Venn2\\IPSC_P3\\B_A_compartment_switch_Venn2.pdf')






#-----------------------------------Venn2_Plot_P20---------------------------------------

data_type = np.dtype({'names':['chr' , 'pos' ,'CCS' , 'NTs' , 'fESC' , 'MEF' , 'IPSC_P3' , 'E14'] , 
                      'formats':['U8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})


union_pc = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NTs = NTs[NTs['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P20 = IPS_P20[IPS_P20['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    for i in range(len(tmp_CCS)):
        union_pc.append((g , i * 200000 , tmp_CCS[i]['pc'] , tmp_NTs[i]['pc'] , tmp_fESC[i]['pc'] , tmp_MEF[i]['pc'] , tmp_IPS_P20[i]['pc'] , tmp_E14[i]['pc']))

union_pc = np.array(union_pc , dtype = data_type)
            




def Common_pc_Venn2(union_pc , lab1 , lab2 , process):
    '''
    '''
    n1 = 0 ; n2 = 0 ; n3 = 0
    
    for i in union_pc:
        if i['CCS'] != 0 and i['NTs'] != 0 and i['fESC'] != 0 and i['MEF'] != 0 and i['IPSC_P3'] != 0 and i['E14'] != 0:
            if process == 'A_B':
                if i['CCS'] > 0 and i['NTs'] < 0 and i['fESC'] < 0 and i['MEF'] > 0 and i['IPSC_P3'] < 0 and i['E14'] < 0:
                    n2 += 1
                elif i['CCS'] > 0 and i['NTs'] < 0 and i['fESC'] < 0 and i['E14'] < 0:
                    n1 += 1
                elif i['MEF'] > 0 and i['IPSC_P3'] < 0 and i['E14'] < 0 and i['fESC'] < 0:
                    n3 += 1
            else:
                if i['CCS'] < 0 and i['NTs'] > 0 and i['fESC'] > 0 and i['MEF'] < 0 and i['IPSC_P3'] > 0 and i['E14'] > 0:
                    n2 += 1
                elif i['CCS'] < 0 and i['NTs'] > 0 and i['fESC'] > 0 and i['E14'] > 0:
                    n1 += 1
                elif i['MEF'] < 0 and i['IPSC_P3'] > 0 and i['E14'] > 0 and i['fESC'] > 0:
                    n3 += 1
                
                

            
    fig = plt.figure(figsize = (12, 12))
    venn2(subsets = (n1, n3, n2), set_labels = (lab1, lab2))  
    
    return fig




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

    
fig1 = Common_pc_Venn2(union_pc , 'A_B_NTs' , 'A_B_IPSC' , 'A_B')
fig2 = Common_pc_Venn2(union_pc , 'B_A_NTs' , 'B_A_IPSC' , 'B_A')

run_Plot(fig1 , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_Venn2\\IPSC_P20\\A_B_compartment_switch_Venn2.pdf')
run_Plot(fig2 , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_Venn2\\IPSC_P20\\B_A_compartment_switch_Venn2.pdf')



#----------------------------------Chromatin_resistant_Venn2_Plot---------------------------------------

data_type = np.dtype({'names':['chr' , 'pos' ,'CCS' , 'NTs' , 'fESC' , 'MEF' , 'IPSC_P3' , 'E14'] , 
                      'formats':['U8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})


union_pc = []
for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NTs = NTs[NTs['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    for i in range(len(tmp_CCS)):
        union_pc.append((g , i * 200000 , tmp_CCS[i]['pc'] , tmp_NTs[i]['pc'] , tmp_fESC[i]['pc'] , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))

union_pc = np.array(union_pc , dtype = data_type)
            




def resistant_pc_Venn2(union_pc , lab1 , lab2 , process):
    '''
    '''
    n1 = 0 ; n2 = 0 ; n3 = 0
    
    for i in union_pc:
        if i['CCS'] != 0 and i['NTs'] != 0 and i['fESC'] != 0 and i['MEF'] != 0 and i['IPSC_P3'] != 0 and i['E14'] != 0:
            if process == 'A_B':
                if i['CCS'] > 0 and i['NTs'] > 0 and i['fESC'] < 0 and i['MEF'] > 0 and i['IPSC_P3'] > 0 and i['E14'] < 0:
                    n2 += 1
                elif i['CCS'] > 0 and i['NTs'] > 0 and i['fESC'] < 0 and i['E14'] < 0:
                    n1 += 1
                elif i['MEF'] > 0 and i['IPSC_P3'] > 0 and i['E14'] < 0 and i['fESC'] < 0:
                    n3 += 1
            else:
                if i['CCS'] < 0 and i['NTs'] < 0 and i['fESC'] > 0 and i['MEF'] < 0 and i['IPSC_P3'] < 0 and i['E14'] > 0:
                    n2 += 1
                elif i['CCS'] < 0 and i['NTs'] < 0 and i['fESC'] > 0 and i['E14'] > 0:
                    n1 += 1
                elif i['MEF'] < 0 and i['IPSC_P3'] < 0 and i['E14'] > 0 and i['fESC'] > 0:
                    n3 += 1
                
                

            
    fig = plt.figure(figsize = (12, 12))
    venn2(subsets = (n1, n3, n2), set_labels = (lab1, lab2))  
    
    return fig




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

    
fig1 = resistant_pc_Venn2(union_pc , 'A_B_NTs' , 'A_B_IPSC' , 'A_B')
fig2 = resistant_pc_Venn2(union_pc , 'B_A_NTs' , 'B_A_IPSC' , 'B_A')

run_Plot(fig1 , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_Venn2\\A_B_compartment_resistant_switch_Venn2.pdf')
run_Plot(fig2 , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_P3_P20_new\\compartment_dynamic_Venn2\\B_A_compartment_resistant_switch_Venn2.pdf')









