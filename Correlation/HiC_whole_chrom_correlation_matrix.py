# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 19:24:48 2018

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys

#from tadlib.calfea import analyze

#--------------------------------------------------------------------------
## Matplotlib Settings
import matplotlib
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF' ,'#CD0000'])
my_cmap.set_bad('#D3D3D3')
                
cells = {'CCS':0 , 'NT5':1 , 'NT6':2 , 'fESC':3}
cell = ['CCS' , 'NT5' , 'NT6' , 'fESC']
chrom=['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19' , 'X']
chrom_1={'1':0 , '2':1 , '3':2 , '4':3 , '5':4 , '6':5 , '7':6 , '8':7 , '9':8 , '10':9 , '11':10 , '12':11 , '13':12 , '14':13 , '15':14 , '16':15 , '17':16 , '18':17 , '19':18 , 'X':19}
HiCFolder = 'D:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_1M'
OutFolder = 'D:\\Workspace_New\\data\\HiC\\Matrix\\matrix_correlation'
size = (12, 12)   
Left = 0.25 ; HB = 0.2 ; width = 0.6 ; HH = 0.6 


matrix_all = {'CCS':[] , 'NT5':[] , 'NT6':[] , 'fESC':[]}

for c in cell:
    wholeFil = c + '_Correct_Merged_Reps_Whole_Genome_Matrix.npz'
    wholeSource = os.path.join(HiCFolder , wholeFil)
    wholeData = np.load(wholeSource)
    matrix = wholeData['Matrix']
    print wholeData['Bins']
    
            
    Bool = np.zeros_like(matrix, dtype = bool)
    for g in chrom:
        start = wholeData['Bins'][()][g][0]
        end = wholeData['Bins'][()][g][1] + 1
        Bool[start:end , start:end] = 1
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if Bool[i][j] == False:
                matrix_all[c].append(matrix[i][j])
            
    

    
cor_matrix = np.zeros((4,4))

for i in cell:
    for j in cell:
        cor_matrix[cells[i] , cells[j]] = np.corrcoef(matrix_all[i] , matrix_all[j])[0][1]
        
o = open(os.path.join(OutFolder , 'whole_chrom_correlation_matrix_1.txt') , 'w')
o.writelines('\t'+'\t'.join(['CCS','NT5','NT6','fESC']) + '\n')
for c in cell:
    o.writelines(c + '\t' + '\t'.join([str(x) for x in cor_matrix[cells[c]]]) + '\n')
o.close()


