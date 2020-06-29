# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 09:18:29 2019

@author: han-luo
"""

import xlrd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def Bar_plot(x , y , colors):
    left, bottom, width, height = 0.1, 0.2, 0.8, 0.6
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.bar(x , y , color = colors)
    xticks = x
    labels = ['CCS' , 'NT5' , 'NT6' , 'fESC']
    for a, b in zip(x, y):
        plt.text(a, b + 0.05, '%.0f' % b, ha='center', va='bottom', fontsize=10)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels,fontsize = 10)
    ax.set_xlabel('Cell type' , fontsize = 20 )
    ax.set_xlim((-1, 4))
    ax.set_ylim((0 , 40000))
    return fig
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

x = [0 , 1 , 2 , 3]
y = [27671 , 34359 , 34338 , 37729]
colors = 'cadetblue'
fig = Bar_plot(x , y , colors)

run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\ATAC_peaknumbers.pdf')