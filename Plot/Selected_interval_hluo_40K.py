# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 11:38:00 2019

@author: Administrator
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
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')
                
def caxis_H(ax):
    """
    Axis Control for HeatMaps.
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = 12, length = 5, pad = 7)
    
def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])
    


#interval = [('2' , 234280000 , 236680000) , ('22' , 33120000 , 35520000) , ('22' , 33320000 , 35720000) , ('22' , 48560000 , 50960000) , 
#            ('22' , 48920000 , 51320000) , ('4' , 7440000 , 9840000) , ('4' , 9360000 , 11760000) , ('4' , 166320000 , 168720000) , 
#            ('7' , 4680000 , 7080000) , ('7' , 4800000 , 7200000) , ('8' , 118920000 , 121320000) , ('X' , 22480000 , 24880000) , 
#            ('X' , 106160000 , 108560000) , ('X' , 108480000 , 110880000)]
        
#interval = [('12' , 5240000 , 7640000) , ('1' , 33000000 , 35400000) , ('10' , 20840000 , 23240000) , ('10' , 116520000 , 118920000) , 
#            ('10' , 116560000 , 118960000) , ('13' , 64120000 , 66520000) , ('12' , 112600000 , 115000000) , ('12' , 113440000 , 115840000) , 
#            ('12' , 113840000 , 116240000) , ('14' , 6440000 , 8840000) , ('14' , 42480000 , 44880000) , ('14' , 42920000 , 45320000) , 
#            ('14' , 43560000 , 45960000) , ('14' , 109760000 , 112160000) , ('17' , 5520000 , 7920000) , ('17' , 17200000 , 19600000) , 
#            ('17' , 18200000 , 20600000) , ('17' , 22080000 , 24480000) , ('17' , 37200000 , 39600000) , ('17' , 83520000 , 85920000) , 
#            ('16' , 43560000 , 45960000) , ('18' , 18400000 , 20800000) , ('18' , 59000000 , 61400000) , ('1' , 170320000 , 172720000) , 
#            ('1' , 172520000 , 174920000) , ('1' , 176600000 , 179000000) , ('3' , 2040000 , 4440000) , ('3' , 93120000 , 95520000) , 
#            ('2' , 19280000 , 21680000) , ('2' , 86520000 , 88920000) , ('5' , 103320000 , 105720000) , ('5' , 103760000 , 106160000) , 
#            ('4' , 144200000 , 146600000) , ('4' , 146640000 , 149400000) , ('7' , 14600000 , 17000000) , ('7' , 40120000 , 42520000) , 
#            ('7' , 40480000 , 42880000) , ('6' , 58320000 , 60720000) , ('6' , 59120000 , 61520000) , ('6' , 59360000 , 61760000) , 
#            ('6' , 59760000 , 62160000) , ('6' , 66360000 , 68760000) , ('6' , 67760000 , 70160000) , ('6' , 69280000 , 71680000) , 
#            ('6' , 129200000 , 131600000) , ('6' , 39000000 , 41400000) , ('9' , 88040000 , 90440000) , ('X' , 7720000 , 10120000) , 
#            ('X' , 36480000 , 38880000) , ('X' , 132840000 , 135240000)]
        

interval = [('6' , 113440000 , 115840000) , ('11' , 58000000 , 60400000) , ('11' , 70120000 , 72520000) , ('13' , 65880000 , 68280000) , 
            ('13' , 112080000 , 114480000) , ('14' , 109760000 , 112160000) , ('16' , 9360000 , 11760000) , ('16' , 89560000 , 91960000) , 
            ('18' , 46480000 , 48880000) , ('1' , 3200000 , 56000000) , ('1' , 7640000 , 10040000) , ('1' , 170320000 , 172720000) , 
            ('1' , 172960000 , 175360000) , ('3' , 30480000 , 32880000) , ('3' , 37680000 , 40080000) , ('2' , 144720000 , 147120000) , 
            ('4' , 111480000 , 113880000) , ('4' , 144200000 , 146600000) , ('7' , 61240000 , 63640000) , ('6' , 131920000 , 134320000) , 
            ('9' , 105280000 , 107680000) , ('9' , 122680000 , 124600000)]
size = (12, 12)
Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
R = 40000

    
f = open('E:\\data\\genome\\mm10.txt' , 'r')
mm ={}
for i in f:
    i = i.split()
    mm[i[0]] = int(i[1])
    

    



tmp = {}




Lib = np.load('D:\\hluo\\Figures\\Matrix\\E7.5_40K_Correct_Merged_Reps_Local_Chromosome_Matrix.npz')

    

pp1 = PdfPages('D:\\hluo\\Figures\\E7.5_40K.pdf')


Lib = Lib['Matrix'][()]

for i in interval:
    print i
    g = i[0]
    startHiC = i[1] // R
    endHiC = i[2] // R
    lib = Lib[g]
    matrix = lib[startHiC:endHiC , startHiC:endHiC]
    nonzero = matrix[np.nonzero(matrix)]
    vmax = np.percentile(nonzero, 95)
        
    #=============HeatMap + colorbar=================================
        
    fig = plt.figure(figsize = size)
    ax1 = fig.add_axes([Left  , HB , width , HH])
    sc = ax1.imshow(matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                    extent = (0, len(matrix), 0, len(matrix)), vmax = vmax, origin = 'lower')
    cxlim = ax1.get_xlim()
    cylim = ax1.get_ylim()
    ## Ticks and Labels
    ticks = list(np.linspace(0 , len(matrix) , 5).astype(float))
    pos = [((startHiC + t) * R) for t in ticks]
    labels = [properU(p) for p in pos]
    ax1.set_xticks(ticks)
    ax1.set_xticklabels(labels)
    ax1.set_yticks(ticks)
    ax1.set_yticklabels(labels, rotation = 'horizontal')
    ax1.set_xlabel('Chr' + i[0] , fontsize = 30 )
        
    ax1.set_xlim(cxlim)
    ax1.set_ylim(cylim)                    
    caxis_H(ax1)
    ## Colorbar
    ax2 = fig.add_axes([Left + 0.6 , HB - 0.06 , 0.1 , 0.035])
    fig.colorbar(sc,cax = ax2, orientation='horizontal' , ticks = [0 , int(vmax/2) , int(vmax)])



    pp1.savefig(fig)
    plt.close(fig)

pp1.close()

