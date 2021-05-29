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
import scipy
from scipy import stats
from scipy.stats import ttest_ind

# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')


def Sig_To_10K(signal , control):
    """
    """
    
    New_Data = {}
    for g in chroms:
        New_Data[g] = {}
        if control == 'input':
            tmp_data = np.array(list(signal.intervals('chr' + g)) , dtype = signal_type)
        else:
            tmp_data = np.array(list(signal.intervals('chr' + g)) , dtype = signal_type)
        max_ = tmp_data['end'].max()
        bin_size = max_ // 10000 + 1
        New_Data[g] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 10000
            New_Data[g][start] += line['value']
    
    return New_Data
    

def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
    
def Box_plot(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(sig['c1'] , 
            positions=[1] , showfliers=False, widths = 0.7 ,
            boxprops={'color': '#A6CEE3','linewidth':0.5},
            medianprops={'color':'#A6CEE3','linewidth':0.5},
            capprops={'color':'#A6CEE3','linewidth':0.5},
            whiskerprops={'color':'#A6CEE3','linewidth':0.5, 'linestyle':'--'})
    ax.boxplot(sig['c2'] ,
            positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': '#1F78B4','linewidth':0.5},
            medianprops={'color':'#1F78B4','linewidth':0.5},
            capprops={'color':'#1F78B4','linewidth':0.5},
            whiskerprops={'color':'#1F78B4','linewidth':0.5, 'linestyle':'--'})
    ax.boxplot(sig['c3'] ,
            positions=[3] , showfliers=False, widths = 0.7, 
            boxprops={'color': '#B2DF8A','linewidth':0.5},
            medianprops={'color':'#B2DF8A','linewidth':0.5},
            capprops={'color':'#B2DF8A','linewidth':0.5},
            whiskerprops={'color':'#B2DF8A','linewidth':0.5, 'linestyle':'--'})
    ax.boxplot(sig['c4'] , 
            positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': '#33A02C','linewidth':0.5},
            medianprops={'color':'#33A02C','linewidth':0.5},
            capprops={'color':'#33A02C','linewidth':0.5},
            whiskerprops={'color':'#33A02C','linewidth':0.5, 'linestyle':'--'})
    ax.boxplot(sig['c5'] , 
            positions=[5] , showfliers=False, widths = 0.7 ,
            boxprops={'color': '#FB9A99','linewidth':0.5},
            medianprops={'color':'#FB9A99','linewidth':0.5},
            capprops={'color':'#FB9A99','linewidth':0.5},
            whiskerprops={'color':'#FB9A99','linewidth':0.5, 'linestyle':'--'})
    ax.boxplot(sig['c6'] , 
            positions=[6] , showfliers=False, widths = 0.7 ,
            boxprops={'color': '#E31A1C','linewidth':0.5},
            medianprops={'color':'#E31A1C','linewidth':0.5},
            capprops={'color':'#E31A1C','linewidth':0.5},
            whiskerprops={'color':'#E31A1C','linewidth':0.5, 'linestyle':'--'})
    ax.boxplot(sig['c7'] , 
            positions=[7] , showfliers=False, widths = 0.7 ,
            boxprops={'color': '#FDBF6F','linewidth':0.5},
            medianprops={'color':'#FDBF6F','linewidth':0.5},
            capprops={'color':'#FDBF6F','linewidth':0.5},
            whiskerprops={'color':'#FDBF6F','linewidth':0.5, 'linestyle':'--'})
    ax.boxplot(sig['c8'] , 
            positions=[8] , showfliers=False, widths = 0.7 ,
            boxprops={'color': '#CAB2D6','linewidth':0.5},
            medianprops={'color':'#CAB2D6','linewidth':0.5},
            capprops={'color':'#CAB2D6','linewidth':0.5},
            whiskerprops={'color':'#CAB2D6','linewidth':0.5, 'linestyle':'--'})
    
#    d1 = scipy.stats.ranksums(data[0] , data[1])[1]
#    d2 = scipy.stats.ranksums(data[0] , data[2])[1]
#    d3 = scipy.stats.ranksums(data[0] , data[3])[1]
#    d4 = scipy.stats.ranksums(data[1] , data[2])[1]
#    d5 = scipy.stats.ranksums(data[1] , data[3])[1]
#    d6 = scipy.stats.ranksums(data[2] , data[3])[1]
    
    ax.set_xticks([1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 , 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23])
    ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F35'  , 'F40' , '' , 'CCs' , 'NT5' , 'NT6' , 'F35' , 'F40'  , '' , 'CCs' , 'NT5' , 'NT6' , 'F35' , 'F40'  , '' , 'CCs' , 'NT5' , 'NT6' , 'F35' , 'F40'] , fontsize = 10)
    ax.set_xlabel('Repr , Partial , Over , Resist')
    ax.set_xlim((0 , 24))
    ax.set_ylim((0 , 40))
    return fig


    
pc_type = np.dtype({'names':['chr' , 'start' , 'end'] , 
                    'formats':['S8' , np.int , np.int]})
signal_type = np.dtype({'names':['start' , 'end' , 'value'] , 
                    'formats':[np.int , np.int , np.float]})

chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000



                    
pc1 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster1_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc2 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster2_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc3 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster3_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc4 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster4_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc5 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster5_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc6 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster6_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc7 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster7_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )
pc8 = np.loadtxt('/public/home/lixinxin/data/BDF1/HiC/Compartment/Compartment_new/Compartment_classify_byself/Compartment_cluster8_pc1.txt' , skiprows = 1 , usecols = (0 , 1, 2) , dtype = pc_type )

# chip1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3_new/mapping/signals/CCs_H3K4me3_R2.bw")
# input1 = pyBigWig.open("/public/home/lixinxin/data/BDF1/Chip/CCS_H3K4me3/mapping/Input.bw")


chip1 = pyBigWig.open("/public/home/lixinxin/data/literature/Gao/Chip/signal_literature/GSM4349664_CC_H3K9me3_rep1.bw")
input1 = pyBigWig.open("/public/home/lixinxin/data/literature/Gao/Chip/signal_literature/GSM4349666_CC_input_rep1.bw")


Chip1 = Sig_To_10K(chip1 , 'chip')
Input1 = Sig_To_10K(input1 , 'input')



pc_data = {'c1' : pc1 , 'c2' : pc2 , 'c3' : pc3, 'c4' : pc4 , 'c5' : pc5 , 'c6' : pc6 , 'c7' : pc7 , 'c8' : pc8 }
sig = {'c1' : [] , 'c2' : [] , 'c3' : [] , 'c4' : [] , 'c5' : [] , 'c6' : [] , 'c7' : [] , 'c8' : []}
cluster = ['c1' , 'c2' , 'c3' , 'c4' , 'c5' , 'c6' , 'c7' , 'c8']



for c in cluster:
    for i in pc_data[c]:
        g = i['chr']
        start = i['start'] // 10000
        end = i['end'] // 10000
        a = Chip1[g][start:end]
        b = Input1[g][start:end]
        for i in range(len(b)):
            if b[i] == 0:
                b[i] = 1
        sig_chip = a / b
        for j in sig_chip:
            sig[c].append(j)
        
        
        
        
pp = PdfPages('/public/home/lixinxin/data/BDF1/Chip/CCS_H3K9me3_Gao/Compartment_cluster_H3K9me3_sig_boxplot_2.pdf')
       
size = (12, 12)
left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
size_axes = [left, bottom, width, height]

color = ['#FB9A99', '#E31A1C', '#FDBF6F', '#CAB2D6', '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',  '#6A3D9A']


fig = plt.figure(figsize = (12, 12))
ax = fig.add_axes(size_axes)
ax.boxplot(sig['c1'] , 
        positions=[1] , showfliers=False, widths = 0.7 ,
        boxprops={'color': color[0],'linewidth':1},
        medianprops={'color':color[0],'linewidth':1},
        capprops={'color':color[0],'linewidth':1},
        whiskerprops={'color':color[0],'linewidth':1, 'linestyle':'--'})
ax.boxplot(sig['c2'] , 
        positions=[2] , showfliers=False, widths = 0.7 ,
        boxprops={'color': color[1],'linewidth':1},
        medianprops={'color':color[1],'linewidth':1},
        capprops={'color':color[1],'linewidth':1},
        whiskerprops={'color':color[1],'linewidth':1, 'linestyle':'--'})
ax.boxplot(sig['c3'] , 
        positions=[3] , showfliers=False, widths = 0.7 ,
        boxprops={'color': color[2],'linewidth':1},
        medianprops={'color':color[2],'linewidth':1},
        capprops={'color':color[2],'linewidth':1},
        whiskerprops={'color':color[2],'linewidth':1, 'linestyle':'--'})
ax.boxplot(sig['c4'] , 
        positions=[4] , showfliers=False, widths = 0.7 ,
        boxprops={'color': color[3],'linewidth':1},
        medianprops={'color':color[3],'linewidth':1},
        capprops={'color':color[3],'linewidth':1},
        whiskerprops={'color':color[3],'linewidth':1, 'linestyle':'--'})
ax.boxplot(sig['c5'] , 
        positions=[6] , showfliers=False, widths = 0.7 ,
        boxprops={'color': color[4],'linewidth':1},
        medianprops={'color':color[4],'linewidth':1},
        capprops={'color':color[4],'linewidth':1},
        whiskerprops={'color':color[4],'linewidth':1, 'linestyle':'--'})
ax.boxplot(sig['c6'] , 
        positions=[7] , showfliers=False, widths = 0.7 ,
        boxprops={'color': color[5],'linewidth':1},
        medianprops={'color':color[5],'linewidth':1},
        capprops={'color':color[5],'linewidth':1},
        whiskerprops={'color':color[5],'linewidth':1, 'linestyle':'--'})
ax.boxplot(sig['c7'] , 
        positions=[8] , showfliers=False, widths = 0.7 ,
        boxprops={'color': color[6],'linewidth':1},
        medianprops={'color':color[6],'linewidth':1},
        capprops={'color':color[6],'linewidth':1},
        whiskerprops={'color':color[6],'linewidth':1, 'linestyle':'--'})
ax.boxplot(sig['c8'] , 
        positions=[9] , showfliers=False, widths = 0.7 ,
        boxprops={'color': color[7],'linewidth':1},
        medianprops={'color':color[7],'linewidth':1},
        capprops={'color':color[7],'linewidth':1},
        whiskerprops={'color':color[7],'linewidth':1, 'linestyle':'--'})

                    
d1 = np.round(ttest_ind(sig['c1'],sig['c8'])[1] , 5)
d2 = np.round(ttest_ind(sig['c2'],sig['c8'])[1] , 5)
d3 = np.round(ttest_ind(sig['c3'],sig['c8'])[1] , 5)
d4 = np.round(ttest_ind(sig['c4'],sig['c8'])[1] , 5)
d5 = np.round(ttest_ind(sig['c5'],sig['c8'])[1] , 5)
d6 = np.round(ttest_ind(sig['c6'],sig['c8'])[1] , 5)
d7 = np.round(ttest_ind(sig['c7'],sig['c8'])[1] , 5)



ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10])
ax.set_xticklabels(['','c1' , 'c2' , 'c3' , 'c4' , '',  'c5' , 'c6' , 'c7' , 'c8',''] ,fontsize = 20)
ax.set_xlabel(['d1:' + str(d1) + ',d2:' + str(d2) + ',d3:' + str(d3) + ',d4:' + str(d4) + ',d5:' + str(d5) + ',d6:' + str(d6) + ',d7:' + str(d7)])
# ax.set_ylim(-5 , 5)
ax.set_ylabel('H3K9me3 signal intensity within compartment')

pp.savefig(fig)
pp.close()             
        
        
 
