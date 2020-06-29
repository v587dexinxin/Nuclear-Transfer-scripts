# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:00:45 2020

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
import pyBigWig
import seaborn as sns
from scipy.interpolate import  interp1d
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


def Sig_To_10K(signal):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        tmp_data = np.array(list(signal.intervals(g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 10000 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 10000
            end = line['end'] // 10000
            for i in range(start,end):
                New_Data[g][i] += line['value']
    
    return New_Data
    
    
    
pc_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                    'formats':['S8' , np.int , np.int]})
signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000


                    
pc1 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster1_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc2 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster2_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc3 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster3_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc4 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster4_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc5 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster5_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc6 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster6_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc7 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster7_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc8 = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/Compartment_cluster8_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )

chip1 = pyBigWig.open("/public/home/xxli/data/literature_data/Gao/Chip/signal/CCS_H3K9me3_chip_50bp.bw")
input1 = pyBigWig.open("/public/home/xxli/data/literature_data/Gao/Chip/signal/CCS_H3K9me3_input_50bp.bw")




Chip1 = Sig_To_10K(chip1)
Input1 = Sig_To_10K(input1)


l = 10000 ; r = 10000
n = (l + res + r) // 10000
pc_data = {'c1' : pc1 , 'c2' : pc2 , 'c3' : pc3, 'c4' : pc4 , 'c5' : pc5 , 'c6' : pc6 , 'c7' : pc7 , 'c8' : pc8 }
sig = {'c1' : np.zeros(n) , 'c2' : np.zeros(n) , 'c3' : np.zeros(n) , 'c4' : np.zeros(n) , 'c5' : np.zeros(n) , 'c6' : np.zeros(n) , 'c7' : np.zeros(n) , 'c8' : np.zeros(n)}
cluster = ['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8']



for c in cluster:
    for i in pc_data[c]:
        g = i['chr']
        start = (i['start'] - l) // 10000 
        end = (i['end'] + r) // 10000
        if end > len(Chip1[g]):
            for j in range(len(Chip1[g]) - start):
                data = Chip1[g][start + j]
                data_input = Input1[g][start + j]
                sums = Chip1[g].sum() 
                sig[c][j] += data / sums
                sig[c][j] -= data_input / sums
#                sig[c][j] -= Input1[g][start + j]
#                sig[c][j] -= Input2[g][start + j]
        else:            
            for j in range(end - start):
                data = Chip1[g][start + j]
                data_input = Input1[g][start + j]
                sums = Chip1[g].sum() 
                sig[c][j] += data / sums
                sig[c][j] -= data_input / sums

    sig[c] = sig[c] / len(pc_data[c])
#                
        
        
 
pp = PdfPages('Compartment_cluster_H3K9me3_ChiP_11.pdf')
       
size = (12, 12)
Left = 0.1 ; HB = 0.15 ; width = 0.6 ; HH = 0.6
color = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#CAB2D6', '#6A3D9A']

fig = plt.figure(figsize = size)
ax = fig.add_axes([Left  , HB , width , HH])
for c in cluster:
    ax.plot(np.array(range(len(sig[c]))), sig[c] , c = color[cluster.index(c)], label = c)
    ax.legend(loc = 'upper right')

pp.savefig(fig)
pp.close()



pp = PdfPages('Compartment_cluster_H3K9me3_ChiP_12.pdf')
       
size = (12, 12)
Left = 0.1 ; HB = 0.15 ; width = 0.6 ; HH = 0.6
color = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#CAB2D6', '#6A3D9A']

fig = plt.figure(figsize = size)
ax = fig.add_axes([Left  , HB , width , HH])
for c in cluster:
    sns.kdeplot(sig[c], kernel='gau' , color = color[cluster.index(c)] , label = c)
    ax.legend(loc = 'upper right')

pp.savefig(fig)
pp.close()



pp = PdfPages('Compartment_cluster_H3K9me3_ChiP_17.pdf')
       
size = (12, 12)
Left = 0.1 ; HB = 0.15 ; width = 0.6 ; HH = 0.6
color = ['#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C', '#FB9A99', '#E31A1C', '#FDBF6F', '#6A3D9A']


fig = plt.figure(figsize = size)
ax = fig.add_axes([Left  , HB , width , HH])
for c in cluster:
    f = interp1d(np.arange(len(sig[c])), sig[c] , kind='cubic')
    x = np.arange(0,len(sig[c]) - 1 , step = 0.01)
    y = f(x)
    ax.plot(x, y , c = color[cluster.index(c)], label = c)
    ax.legend(loc = 'upper right')

pp.savefig(fig)
pp.close()
