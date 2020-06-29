# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 17:04:16 2020

@author: han-luo
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
            


data1 = pd.read_excel('D:\\ntESC_3Dreprogramming\\Paper\\12345678 GO.xlsx' , index_col=1 , sheetname=0)
data2 = pd.read_excel('D:\\ntESC_3Dreprogramming\\Paper\\12345678 GO.xlsx' , index_col=1 , sheetname=1)
data3 = pd.read_excel('D:\\ntESC_3Dreprogramming\\Paper\\12345678 GO.xlsx' , index_col=1 , sheetname=2)
data4 = pd.read_excel('D:\\ntESC_3Dreprogramming\\Paper\\12345678 GO.xlsx' , index_col=1 , sheetname=3)
data5 = pd.read_excel('D:\\ntESC_3Dreprogramming\\Paper\\12345678 GO.xlsx' , index_col=1 , sheetname=4)
data6 = pd.read_excel('D:\\ntESC_3Dreprogramming\\Paper\\12345678 GO.xlsx' , index_col=1 , sheetname=5)
data7 = pd.read_excel('D:\\ntESC_3Dreprogramming\\Paper\\12345678 GO.xlsx' , index_col=1 , sheetname=6)
data8 = pd.read_excel('D:\\ntESC_3Dreprogramming\\Paper\\12345678 GO.xlsx' , index_col=1 , sheetname=7)


a = data1.iloc[[3 , 2 , 1 , 0],:].PValue
b = data2.iloc[[9 , 8 , 5 , 2 , 0],:].PValue
c = data3.PValue
d = data4.iloc[[9 , 6 , 2 , 0],:].PValue
e = data5.iloc[[10 , 5 , 3 , 0],:].PValue
f = data6.iloc[[1, 0],:].PValue
g = data7.iloc[[9 , 3 , 0],:].PValue
h = data8.iloc[[4 , 0],:].PValue


y1 = -np.log10(np.array(h))
y2 = -np.log10(np.array(g))
y3 = -np.log10(np.array(f))
y4 = -np.log10(np.array(e))

left, bottom, width, height = 0.45, 0.1, 0.5, 0.8
size_axes = [left, bottom, width, height]
fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
y = list(y1) + [0] + list(y2) + [0] + list(y3) + [0] + list(y4) 
colors = ['orange' for i in range(len(y1))] + ['orange'] + ['skyblue' for i in range(len(y2))] + ['skyblue'] + ['greenyellow' for i in range(len(y3))] + ['greenyellow'] + ['magenta' for i in range(len(y4))]
ax.barh(range(len(y1) + len(y2) + len(y3) + len(y4) + 3) , y , color = colors)

yticks = range(len(y1) + len(y2) + len(y3) + len(y4) + 3)
labels = [i.split('~')[1] for i in h.index] + [''] + [i.split('~')[1] for i in g.index] + [''] + [i.split('~')[1] for i in f.index] + [''] + [i.split('~')[1] for i in e.index] 
ax.set_yticks(yticks)
ax.set_yticklabels(labels,fontsize = 10)
ax.set_xlabel('-log10(p value)' , fontsize = 20 )

run_Plot(fig ,'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S4_figs\\Compartment_B_A_GO.pdf' )
