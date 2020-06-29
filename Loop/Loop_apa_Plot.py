# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 10:09:20 2019

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
from scipy.special import ndtr
import math
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


    
def apa_submatrix(M, pos, w=5):
    
    Len = M.shape[0]

    apa = []
    for i, j in pos:
        if (i-w>=0) and (i+w+1<=Len) and (j-w>=0) and (j+w+1<=Len):
            tmp = M[i-w:i+w+1, j-w:j+w+1]
            if tmp.mean()==0:
                continue
            mask = np.isnan(tmp)
            if mask.sum() > 0:
                continue
            tmp = tmp / tmp.mean()
            apa.append(tmp)
    
    return apa

def apa_analysis(apa, w=5, cw=3):
    
    avg = apa.mean(axis=0)
    lowerpart = avg[-cw:,:cw]
    upperpart = avg[:cw,-cw:]
    maxi = upperpart.mean() * 5
    ## APA score
    score = avg[w,w] / lowerpart.mean()
    ## z-score
    z = (avg[w,w] - lowerpart.mean()) / lowerpart.std()
    p = 1 - ndtr(z)
    
    return avg, score, z, p, maxi
    
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

chros = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19']
res = 20000
initial_distance = 2 * math.sqrt(2)

data_type = np.dtype({'names':['chr' , 'start' , 'end'] ,
                      'formats':['S8' , np.int , np.int]})

cell = 'fESC'
#Matrix_data
Matrix = np.load('/public/home/xxli/data/BDF1_New/HiC/HapHiC_workspace/' + cell + '_workspace/' + cell + '_20K/Correct_Merged_Reps_Local_Chromosome_Matrix.npz',allow_pickle=True)
Matrix = Matrix['Matrix'][()]

#Loop_data
P_Loops = np.loadtxt('/public/home/xxli/data/BDF1_New/HiC/Loop/Raw_20K_0.05/Cluster_' + cell + '_loops20K_3.txt' , 
                     skiprows = 1 , dtype = data_type , usecols = (0 , 1 , 2))
                     

apa = []                    
for c in chros:
    peaks = P_Loops[P_Loops['chr'] == c]
    pos = []
    M = Matrix[c]
    for p in peaks:
        x, y = p[1], p[2]
        if abs(y-x) < 15 * res:
            continue
        s_l = range(p[1]//res, int(np.ceil((p[1]+20000)/float(res))))
        e_l = range(p[2]//res, int(np.ceil((p[2]+20000)/float(res))))
        si, ei = None, None
        for st in s_l:
            for et in e_l:
                if (st < M.shape[0]) and (et < M.shape[0]):
                    if si is None:
                        si, ei = st, et
                    else:
                        if M[st,et]>M[si,ei]:
                            si,ei = st,et
        if not si is None:
            if si < ei:
                pos.append((si,ei))
            else:
                pos.append((ei,si))
    tmp = apa_submatrix(M, pos)
    apa.extend(tmp)
apa = np.r_[apa]
avg,score,z,p,maxi = apa_analysis(apa)


Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6
fig = plt.figure(figsize = (12,12))
plt.tick_params(axis='both', bottom=False, top=False, left=False, right=False,
                labelbottom=False, labeltop=False, labelleft=False, labelright=False)
ax = fig.add_axes([Left  , HB , width , HH])
sc = ax.imshow(avg, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
               extent = (0, len(avg), 0, len(avg)), vmin = 0 , vmax = 3, origin = 'lower')

ax.set_xticks([0 , 5 , 10])
ax.set_xticklabels(['-100K' , 'Loop'  , '100K'])
ax.set_yticks([0 , 5 , 10])
ax.set_yticklabels(['-100K' , 'Loop'  , '100K'])
ax.set_title('APA score = {0:.3g}, p-value = {1:.3g}'.format(score, p))
ax = fig.add_axes([Left + width + 0.03 , HB , 0.035 , 0.1])
cbar = fig.colorbar(sc,cax = ax, orientation='vertical')
cbar.set_ticks([0 , 3])

run_Plot(fig , cell + '_Loops_apa_20K.pdf')