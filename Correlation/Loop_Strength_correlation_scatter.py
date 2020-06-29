# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 15:55:03 2018

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import os
import sys


c = 'fESC'
data_type = np.dtype({'names':['chr' , 'R1' , 'R2'],
                       'formats':['S4' , np.float , np.float]})
f = np.loadtxt('D:\\Workspace_New\\data\\HiC\\Loop\\' + c + '_R1_R2_point.txt', 
               dtype = data_type , skiprows = 1, usecols = (0 , 3 , 4))
x = []
y = []

for i in f:
    if (i['R1'] > 1) and (i['R2'] > 1) and (i['R1'] < 50) and (i['R2'] < 50):
        x.append(i['R1'])
        y.append(i['R2'])
        


size = (12, 12)
Left = 0.2 ; HB = 0.2 ; width = 0.6 ; HH = 0.6

        
pp = PdfPages('D:\\Workspace_New\\Plot\\Correlation\\replicated_scatter\\Loop_strength\\' + c + 'loop_strength_R1_R2.pdf')
cor = round(np.corrcoef(x , y)[0][1],5)
fig = plt.figure(figsize = (10, 10))
ax = fig.add_axes([Left , HB , width, HH])
ax.scatter(x , y , alpha = 0.6 , c = 'red')
ax.set_xlabel(c + '_R1' , size = 50)
ax.set_ylabel(c + '_R2' , size = 50)
ax.set_xlim(0,50)
ax.set_ylim(0,50)
ax.text(5 , 40 , 'R = ' + str(cor) , size = 50 )
#    ax.plot([-0.13 , 0.13] , [0 , 0] , ls = '--' , c = 'black' , lw = 1.0 )
#    ax.plot([0 , 0] , [-0.13 , 0.13] , ls = '--' , c = 'black' , lw = 1.0 )
pp.savefig(fig)
pp.close()       