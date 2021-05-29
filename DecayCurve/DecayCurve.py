# -*- coding: utf-8 -*-
"""
Created on Mon Nov 04 10:31:56 2019

@author: han-luo
"""

from __future__ import division
import numpy as np
import sys
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages




def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
    # load HiC data
CCS = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Cooler_to_Matrix\\sparse\\spase_40K\\Corrected_CCS_40K.npz'
CCS = np.load(CCS)
NT5 = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Cooler_to_Matrix\\sparse\\spase_40K\\Corrected_NT5_40K.npz'
NT5 = np.load(NT5)
NT6 = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Cooler_to_Matrix\\sparse\\spase_40K\\Corrected_NT6_40K.npz'
NT6 = np.load(NT6)
fESC = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Cooler_to_Matrix\\sparse\\spase_40K\\Corrected_fESC_40K.npz'
fESC = np.load(fESC)
F35 = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Cooler_to_Matrix\\sparse\\spase_40K\\Corrected_F35_40K.npz'
F35 = np.load(F35)
HiCLib = {'CCS':CCS , 'NT5':NT5 , 'NT6':NT6 , 'fESC':fESC , 'F35':F35}
cell = ['CCS' , 'NT5' , 'NT6' ,'fESC' , 'F35']

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
        print ("chromosome %s is done" % g)
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
ax.plot(range(2,250), data[3][1:] / sum(data[3][1:]), label = 'F40')
ax.plot(range(2,250), data[4][1:] / sum(data[4][1:]), label = 'F35')

#ax.set_xticks(chrlen)
#ax.set_xticklabels(chromosome)
ax.set_xlabel('Distance(40Kb)')
#ax.set_xlim([0,50])
#ax.set_yticks([0.5, 0.75, 0.9])
#ax.set_yticklabels(['0.5', '0.75', '0.9'])
ax.set_ylabel('Contact frequency')
ax.legend()

plt.xscale('log')
plt.yscale('log')

run_Plot(fig , 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Decaycurve\\DecayCurve_all_chrom_normalization.pdf')


