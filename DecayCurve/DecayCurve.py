# -*- coding: utf-8 -*-
"""
Created on Mon Nov 04 10:31:56 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
import sys, cPickle
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
    # load HiC data
CCS = 'H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\sparse\\Corrected_CCS_40K.npz'
CCS = np.load(CCS)
NT5 = 'H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\sparse\\Corrected_NT5_40K.npz'
NT5 = np.load(NT5)
NT6 = 'H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\sparse\\Corrected_NT6_40K.npz'
NT6 = np.load(NT6)
fESC = 'H:\\Workspace_New\\data\\HiC\\Matrix\\Matrix_40K\\sparse\\Corrected_fESC_40K.npz'
fESC = np.load(fESC)
HiCLib = {'CCS':CCS , 'NT5':NT5 , 'NT6':NT6 , 'fESC':fESC}
cell = ['CCS' , 'NT5' , 'NT6' ,'fESC']

    #resolution
Res = 40000
    #chromosome
chrmo = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
#chrmo = [i for i in HiClib.keys() if i != 'resolution' and i != 'genomeInformation']
def distance_contact(HiClib):
    distance_contact = {}
    #max distance 10Mb
    max_distance = 10000000
    #distance array step : Res
    distance = np.arange(1,max_distance//Res,1)

    for i in distance:
        distance_contact[i] = []

    for g in chrmo:
        # contact Matrix
        Matrix = HiClib[g] 
        for i in distance:
            mask = (Matrix['bin2'] - Matrix['bin1']) == i
            IF = Matrix[mask]['IF'].mean()
            distance_contact[i] += [IF]
        print "chromosome %s is done" % g
    for dis,value in distance_contact.items():
        distance_contact[dis] = np.array(value).mean()
    return distance_contact
    
data = []
for c in cell:    
    distance = distance_contact(HiCLib[c])
    y = [distance[i] for i in range(1 , 250)]
    data.append(y)
    
    


fig = plt.figure(figsize = (12, 6))
left, width, bottom, height = 0.2, 0.60, 0.2, 0.60
size_axes = [left, bottom, width, height]
ax = fig.add_axes(size_axes)
ax.plot(range(2,250), data[0][1:] / sum(data[0][1:]), label = 'CCS')
ax.plot(range(2,250), data[1][1:] / sum(data[1][1:]), label = 'NT5')
ax.plot(range(2,250), data[2][1:] / sum(data[2][1:]), label = 'NT6')
ax.plot(range(2,250), data[3][1:] / sum(data[3][1:]), label = 'fESC')

#ax.set_xticks(chrlen)
#ax.set_xticklabels(chromosome)
ax.set_xlabel('Distance(40Kb)')
#ax.set_xlim([0,50])
#ax.set_yticks([0.5, 0.75, 0.9])
#ax.set_yticklabels(['0.5', '0.75', '0.9'])
ax.set_ylabel('Contact frequency')
ax.legend()

#plt.xscale('log')
plt.yscale('log')

run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\Plot\\S4_figs\\DecayCurve_all_chrom_normalization.pdf')


