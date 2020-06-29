# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 20:55:56 2020

@author: han-luo
"""


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
    
Left = 0.15 ; HB = 0.15 ; width = 0.7 ; HH = 0.7
fig = plt.figure(figsize = (14, 8))
ax = fig.add_axes([Left  , HB , width , HH])
ax.scatter(range(1,9) , [0.1448087 , 0.1 , 0.2 , 0.0947368 , 0.39666239 , 0.45 , 0.44 , 0.73584906 ] ,s = 50)
ax.plot(range(1,9) , [0.1448087 , 0.1 , 0.2 , 0.0947368 , 0.39666239 , 0.45 , 0.44 , 0.73584906 ])
ax.set_xlim((0,9))
ax.set_ylim((0,1))
ax.set_xticks([1,2,3,4,5,6,7,8])
ax.set_xticklabels(['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8'] ,fontsize = 20)
ax.set_ylabel('Percentage' , fontsize = 30)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig6_figs\\H3K9me3_marked_cluster_percentage.pdf')