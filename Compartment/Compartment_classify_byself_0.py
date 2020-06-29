# -*- coding: utf-8 -*-
"""
Created on Sat Nov 09 15:36:49 2019

@author: han-luo
"""

from __future__ import division
from scipy import sparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import cPickle
import sys
import os

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['darkblue','darkorange'])
my_cmap.set_bad('#2672a1')


pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['S8' , np.float]})
                    
                    
def Classify3_compartment(comp1 , comp2 , comp3):
    classify = {'AAA':[] , 'AAB':[] , 'ABA':[] , 'ABB':[] , 'BBB':[] , 'BBA':[] , 'BAA':[] , 'BAB':[]}
    for g in chrom:
        tmp_comp1 = comp1[comp1['chr'] == g]
        tmp_comp2 = comp2[comp2['chr'] == g]
        tmp_comp3= comp3[comp3['chr'] == g]
        for i in range(len(tmp_comp1)):
            score1 = tmp_comp1[i]['pc']
            score2 = tmp_comp2[i]['pc']
            score3 = tmp_comp3[i]['pc']
            if (score1 > 0) and (score2 > 0) and (score3 > 0):
                classify['AAA'].append((score1 , score2 , score3))
            elif (score1 > 0) and (score2 > 0) and (score3 < 0):
                classify['AAB'].append((score1 , score2 , score3))
            elif (score1 > 0) and (score2 < 0) and (score3 > 0):
                classify['ABA'].append((score1 , score2 , score3))
            elif (score1 > 0) and (score2 < 0) and (score3 < 0):
                classify['ABB'].append((score1 , score2 , score3))
            elif (score1 < 0) and (score2 < 0) and (score3 < 0):
                classify['BBB'].append((score1 , score2 , score3))
            elif (score1 < 0) and (score2 < 0) and (score3 > 0):
                classify['BBA'].append((score1 , score2 , score3))
            elif (score1 < 0) and (score2 > 0) and (score3 > 0):
                classify['BAA'].append((score1 , score2 , score3))
            elif (score1 < 0) and (score2 > 0) and (score3 < 0):
                classify['BAB'].append((score1 , score2 , score3))
    return classify
    
    
CCS = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\CCS_compartment_200K.txt' , dtype = pc_type)
CCS = CCS[CCS['chr'] != 'X']
NT5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT5_compartment_200K.txt' , dtype = pc_type)
NT5 = NT5[NT5['chr'] != 'X']
NT6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\NT6_compartment_200K.txt' , dtype = pc_type)
NT6 = NT6[NT6['chr'] != 'X']
fESC = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Compartment\\compartment_new\\fESC_compartment_200K.txt' , dtype = pc_type)
fESC = fESC[fESC['chr'] != 'X']

MEF = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_Compartment_200K.txt' , dtype = pc_type)
MEF = MEF[MEF['chr'] != 'X']
IPS_P3 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_Compartment_200K.txt' , dtype = pc_type)
IPS_P3 = IPS_P3[IPS_P3['chr'] != 'X']
E14 = np.loadtxt('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_Compartment_200K.txt' , dtype = pc_type)
E14 = E14[E14['chr'] != 'X']


chrom = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']

##mean_NT5_NT6

NT = []
for g in chrom:
    tmp_NT5 = NT5[NT5['chr'] == g]
    tmp_NT6 = NT6[NT6['chr'] == g]
    for i in range(len(tmp_NT5)):
        score1 = tmp_NT5[i]['pc']
        score2 = tmp_NT6[i]['pc']
        score = (score1 + score2) / 2
        NT.append((g , score))
    
NT = np.array(NT , dtype = pc_type)

classify1 = Classify3_compartment(CCS , NT , fESC)
classify2 = Classify3_compartment(MEF , IPS_P3 , E14)

    
matrix = np.array([[1,1,1] , [1,1,-1] , [1,-1,-1] , [-1,-1,-1] , [-1,1,1] , [-1,1,-1]])
size = (12, 12)
Left = 0.35 ; HB = 0.15 ; width = 0.35 ; HH = 0.7

fig = plt.figure(figsize = size)
ax1 = fig.add_axes([Left  , HB , width , HH])
sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', vmax = 1 , vmin = -1 , extent = (0, 3, len(matrix) , 0))    
    
            
for i in range(1 , 7):
    ax1.plot([0 , 3] , [i , i] , color="black")
for i in range(1 , 3):
    ax1.plot([i , i] , [0 , 6] , color="black")

                