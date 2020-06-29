# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 19:38:40 2018

@author: xxli
"""

import numpy as np  
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


Left = 0.3 ; HB = 0.15 ; width = 0.4 ; HH = 0.7

my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])

def DrawBubble(read_name):#气泡图
    sns.set(style = "whitegrid")#设置样式
    
    x = range(2) * 20
    y = []
    for i in range(20):
        y.extend([i,i])
    
    z = []
    for i in range(len(x)):
        z.append(p_matrix[y[i]][x[i]])#用来调整各个点的大小s
    fpkm = []
    for i in range(len(x)):
        fpkm.append(g_matrix[y[i]][x[i]])#用来调整各个点的大小s
    cm = plt.cm.get_cmap('Reds')
    fig = plt.figure(figsize = (10 , 12))
    ax = fig.add_axes([Left  , HB , width , HH])
    bubble = ax.scatter(x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5, alpha = 0.8)
#    ax.scatter([x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5, alpha = 0.8)
    ax.grid()
    
    ax.set_xticks([-1, 0 , 1 , 2 ])
    ax.set_xticklabels(['' , 'CCS' , 'NT_fESC' , ''] , fontsize = 15)
    ax.set_yticks(range(20))
    ax.set_yticklabels(gene_name , fontsize = 20)
    ax.set_xlabel('classifies', fontsize = 25 , labelpad = 40)#X轴标签
    ax.set_ylabel('TFs', fontsize = 25)#Y轴标签
    
    ax2 = fig.add_axes([Left + 0.34 , HB - 0.05 , 0.06 , 0.015])
    fig.colorbar(bubble,cax = ax2, orientation='horizontal' , ticks = [np.round(min(fpkm),2) , int(max(fpkm)/2) , int(max(fpkm))])
    
    
    plt.savefig('D:\\Workspace_New\\Plot\\Figure3\\motif_enrichment.png')
if __name__=='__main__':
    DrawBubble("PeopleNumber.csv")#气泡图


def Bar_Plot():
    Left = 0.1 ; HB = 0.15 ; width = 0.8 ; HH = 0.7
    bar_width = 1
    opacity = 0.8
    x0 = np.arange(1 , 27  , 5)
    y0 = np.array([g1['CCS'][0] ,g2['CCS'][0] ,g3['CCS'][0] ,g4['CCS'][0] ,g5['CCS'][0] ,g6['CCS'][0]])
    x1 = x0 + 1
    y1 = np.array([g1['NT5'][0] ,g2['NT5'][0] ,g3['NT5'][0] ,g4['NT5'][0] ,g5['NT5'][0] ,g6['NT5'][0]])
    x2 = x1 + 1
    y2 = np.array([g1['NT6'][0] ,g2['NT6'][0] ,g3['NT6'][0] ,g4['NT6'][0] ,g5['NT6'][0] ,g6['NT6'][0]])
    x3 = x2 + 1
    y3 = np.array([g1['fESC'][0] ,g2['fESC'][0] ,g3['fESC'][0] ,g4['fESC'][0] ,g5['fESC'][0] ,g6['fESC'][0]])
    fig = plt.figure(figsize = (12, 10))
    ax = fig.add_axes([Left  , HB , width , HH])
    ax.bar(x0, y0, bar_width,alpha=opacity, color='b',label=    'CCS')
    ax.bar(x1, y1, bar_width,alpha=opacity, color='r',label=    'NT5')
    ax.bar(x2, y2, bar_width,alpha=opacity, color='g',label=    'NT6')
    ax.bar(x3, y3, bar_width,alpha=opacity, color='y',label=    'fESC')
    ax.set_xticks(x1)
    ax.set_xticklabels(['Atm' , 'Atr' , 'Trp53' , 'Rad50' , 'Parp1' , 'Brca1'] , fontsize = 20)
    ax.set_xlim((0 , 30))
    ax.set_yticks([0 , 50 , 100 , 150 , 200 , 250])
    ax.set_yticklabels([0 , 50 , 100 , 150 , 200 , 250] , fontsize = 15)
    ax.set_ylabel('FPKM' , fontsize = 25)
    ax.legend()


def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    


Left = 0.3 ; HB = 0.15 ; width = 0.4 ; HH = 0.7
cm = plt.cm.get_cmap('Reds')
fig = plt.figure(figsize = (10 , 12))
ax = fig.add_axes([Left  , HB , width , HH])
bubble = ax.scatter(x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5 , alpha = 0.8)
#ax.scatter([x, y , s = np.array(z) / 2 + 5, c = np.array(fpkm) + 0.5, cmap = cm, linewidth = 0.5, alpha = 0.8)
ax.grid()
    
ax.set_xticks([-1, 0 , 1 , 2 ])
ax.set_xticklabels(['' , 'CCS' , 'NT_fESC' , ''] , fontsize = 15)
ax.set_yticks(range(20))
ax.set_yticklabels(gene_name , fontsize = 20)
ax.set_xlabel('classifies', fontsize = 25 , labelpad = 40)#X轴标签
ax.set_ylabel('TFs', fontsize = 25)#Y轴标签

ax2 = fig.add_axes([Left + 0.34 , HB - 0.05 , 0.06 , 0.015])
fig.colorbar(bubble,cax = ax2, orientation='horizontal' , ticks = [np.round(min(fpkm),1) + 0.6, int(max(fpkm)/2) , int(max(fpkm))])

    
    
ax = fig.add_axes([Left + width + 0.05  , HB , 0.2 , 0.3])
ax.scatter([1,1,1,1,1] , [1,2,3,4,5] ,s = [5 , 130 , 255 , 380 , 505] , marker = 'o' , c = 'white' ,linewidth = 0.5,edgecolors = 'black')
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\promoter_classify_motif_bobble.pdf')