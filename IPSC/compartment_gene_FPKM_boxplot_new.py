# -*- coding: utf-8 -*-
"""
Created on Wed May 06 14:24:57 2020

@author: han-luo
"""

from __future__ import division
import math
import numpy as np
import csv , copy
import xlrd
from itertools import islice
from sklearn.cluster import KMeans
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster import hierarchy
from scipy import cluster 
import copy
import scipy
from scipy import stats
import scipy.cluster.hierarchy as sch
from itertools import islice  
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


def get_PC(PC_fil):
    '''
    '''
    pc_type_1 = np.dtype({'names':['chr' , 'pos' , 'pc'] , 
                        'formats':['U8' , np.int , np.float]})
    
    
    chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
    PC_new = []
    PC = np.loadtxt(PC_fil , dtype = pc_type)
    PC = PC[PC['chr'] != 'X']
    for g in chroms:
        tmp = PC[PC['chr'] == g]
        for i in range(len(tmp)):
            PC_new.append((g , i * 200000 , tmp[i]['pc']))
    PC_new = np.array(PC_new , dtype = pc_type_1)
    
    return PC_new


def get_union_gene_sites(Fil):
    '''
    '''
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F40_FPKM'],
                 'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float , np.float]})
         
                 
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if (i['CCS_FPKM'] > 0.1) or (i['NT5_FPKM'] > 0.1) or (i['NT6_FPKM'] > 0.1) or (i['F35_FPKM'] > 0.1) or (i['F40_FPKM'] > 0.1):
            gene_site = int((i['start'] + i['end']) / 2)
            union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS_FPKM'] ,  i['NT5_FPKM'] , i['NT6_FPKM'] , i['F35_FPKM'] , i['F40_FPKM']))
    union = np.array(union , dtype = union_type_1)
    return union

def get_union_gene_all_site_IPSC(Fil):
    union_type = ({'names':['gene_id' , 'gene_name','chr','strand','start','end','CCS_1','CCS_2','CCS_3','NT5_1','NT5_2','NT5_3','NT5_4','NT6_1','NT6_2','NT6_3','F35_1','F35_2','F35_3','F40_1','F40_2','F40_3','MEF_1','MEF_2','IPSC_1','IPSC_2','E14_1','E14_2'],
                 'formats':['U64' , 'U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float , np.float , np.float, np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_id','gene_name','chr','pos','CCS','NT5','NT6','F35','F40','MEF','IPSC','E14'],
                 'formats':['U64' , 'U64' , 'U8' , np.int , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        gene_site = int((i['start'] + i['end']) / 2)
        ccs = (i['CCS_1'] + i['CCS_2'] +i['CCS_3']) / 3
        nt5 = (i['NT5_1'] + i['NT5_2'] +i['NT5_3']) / 3
        nt6 = (i['NT6_1'] + i['NT6_2'] +i['NT6_3']) / 3
        f35 = (i['F35_1'] + i['F35_2'] +i['F35_3']) / 3
        f40 = (i['F40_1'] + i['F40_2'] +i['F40_3']) / 3
        mef = (i['MEF_1'] + i['MEF_2']) / 2
        ipsc = (i['IPSC_1'] + i['IPSC_2']) / 2
        e14 = (i['E14_1'] + i['E14_2']) / 2
        if (ccs > 0.1) or (nt5 > 0.1) or (nt6 > 0.1) or (f35 > 0.1) or (f40 > 0.1) or (mef > 0.1) or (ipsc > 0.1) or (e14 > 0.1):
            tmp = (i['gene_id'] , i['gene_name'] , i['chr'] , gene_site , ccs , nt5 , nt6 , f35 , f40 , mef , ipsc , e14)
            union.append(tmp)
    union = np.array(union , dtype = union_type_1)
    return union

def get_union_gene_all_site_IPSC(Fil):
    union_type = ({'names':['gene_name','chr','strand' , 'start','end','CCS','NT5','NT6','F35','F40','MEF','IPSC','E14'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_name','chr','strand' , 'pos','CCS','NT5','NT6','F35','F40','MEF','IPSC','E14'],
                     'formats':['U64' , 'U8' , 'U8' , np.int , np.float , np.float , np.float , np.float, np.float , np.float , np.float , np.float]})    
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if (i['CCS'] > 0.1) or (i['NT5'] > 0.1) or (i['NT6'] > 0.1) or (i['F35'] > 0.1) or (i['F40'] > 0.1) or (i['MEF'] > 0.1) or (i['IPSC'] > 0.1) or (i['E14'] > 0.1):
            union.append((i['gene_name'] , i['chr'] , i['strand'] , (i['start'] + i['end']) / 2 , i['CCS'] , i['NT5'] , i['NT6'] , i['F35'] , i['F40'] , i['MEF'] , i['IPSC'] , i['E14']))
    union = np.array(union , dtype = union_type_1)
    return union




def get_raw_genes_new(Fil):
    '''
    '''
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT2' , 'NT3' , 'NT4' , 'NT5' , 'NT6' , 'F35' , 'F37' , 'F40' , 'F41'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        if (i['CCS'] > 0.1) or (i['NT5'] > 0.1) or (i['NT6'] > 0.1) or (i['F35'] > 0.1) or (i['F40'] > 0.1):
            union.append(i)
    union = np.array(union , dtype = union_type)

    return union
    


def Load_gtf(gtfil):
    '''
    '''
    gtf_type = np.dtype({'names':['gene_id' , 'gene_name' , 'chr' , 'strand' , 'start' , 'end'],
                     'formats':['U64' , 'U64' , 'U8' , 'U4' , np.int , np.int]})
    gtf = open(gtfil , 'r')
    gtf_1 = []
    for i in islice(gtf , 5 , None):
        a = i.strip().split()
        if a[2] == 'gene':
            gene_id = i.strip().split('\"')[1]
            gene_name = i.strip().split('\"')[5]
            chro = a[0]
            strand = a[6]
            start = a[3]
            end = a[4]
            gtf_1.append((gene_id , gene_name , chro , strand , start , end))
    gtf = np.array(gtf_1 , dtype = gtf_type)
    return gtf

               
    
def Sort(a , s_1):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1]))
    return a



 
# def Box_plot(data):                
#     left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
#     size_axes = [left, bottom, width, height]
#     fig = plt.figure(figsize = (12, 12))
#     ax = fig.add_axes(size_axes)
#     ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 ,
#             boxprops={'color': 'seagreen','linewidth':2},
#             medianprops={'color':'seagreen','linewidth':2},
#             capprops={'color':'seagreen','linewidth':2},
#             whiskerprops={'color':'seagreen','linewidth':2, 'linestyle':'--'})
#     ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
#             boxprops={'color': 'chocolate','linewidth':2},
#             medianprops={'color':'chocolate','linewidth':2},
#             capprops={'color':'chocolate','linewidth':2},
#             whiskerprops={'color':'chocolate','linewidth':2, 'linestyle':'--'})
#     ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7, 
#             boxprops={'color': 'slateblue','linewidth':2},
#             medianprops={'color':'slateblue','linewidth':2},
#             capprops={'color':'slateblue','linewidth':2},
#             whiskerprops={'color':'slateblue','linewidth':2, 'linestyle':'--'})
#     ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
#             boxprops={'color': 'deeppink','linewidth':2},
#             medianprops={'color':'deeppink','linewidth':2},
#             capprops={'color':'deeppink','linewidth':2},
#             whiskerprops={'color':'deeppink','linewidth':2, 'linestyle':'--'})
#     ax.boxplot(data[4] , positions=[5] , showfliers=False, widths = 0.7 ,
#             boxprops={'color': 'darkmagenta','linewidth':2},
#             medianprops={'color':'darkmagenta','linewidth':2},
#             capprops={'color':'darkmagenta','linewidth':2},
#             whiskerprops={'color':'darkmagenta','linewidth':2, 'linestyle':'--'})

# #    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
#     d1 = scipy.stats.ranksums(data[0] , data[1])[1]
#     d2 = scipy.stats.ranksums(data[0] , data[2])[1]
#     d3 = scipy.stats.ranksums(data[0] , data[3])[1]
#     d4 = scipy.stats.ranksums(data[0] , data[4])[1]
#     d5 = scipy.stats.ranksums(data[1] , data[2])[1]
#     d6 = scipy.stats.ranksums(data[1] , data[3])[1]
#     d7 = scipy.stats.ranksums(data[1] , data[4])[1]
#     d8 = scipy.stats.ranksums(data[2] , data[3])[1]
#     d9 = scipy.stats.ranksums(data[2] , data[4])[1]
#     d10 = scipy.stats.ranksums(data[3] , data[4])[1]
    
    
#     ax.set_xticks([1 , 2 , 3 , 4])
#     ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F35' , 'F40' ] , fontsize = 28)
#     ax.set_xlabel('CCS_NT5:' + str(d1) + ',CCS_NT6:' + str(d2) + ',CCS_F35:' + str(d3) + ',CCS_F40:' + str(d4) + 'NT5_NT6:' + str(d5) + '\n' +
#                   ',NT5_F35:' + str(d6) + ',NT5_F40:' + str(d7) + ',NT6_F35:' + str(d8) + ',NT6_F40:' + str(d9) + ',F35_F40:' + str(d10))
#     ax.set_xlim((0.5 , 5.5))
#     ax.set_ylim((0 , 25))
#     return fig



def Box_plot(data , label , vmin , vmax):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'seagreen','linewidth':1},
            medianprops={'color':'seagreen','linewidth':1},
            capprops={'color':'seagreen','linewidth':1},
            whiskerprops={'color':'seagreen','linewidth':1, 'linestyle':'--'})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'chocolate','linewidth':1},
            medianprops={'color':'chocolate','linewidth':1},
            capprops={'color':'chocolate','linewidth':1},
            whiskerprops={'color':'chocolate','linewidth':1, 'linestyle':'--'})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7, 
            boxprops={'color': 'slateblue','linewidth':1},
            medianprops={'color':'slateblue','linewidth':1},
            capprops={'color':'slateblue','linewidth':1},
            whiskerprops={'color':'slateblue','linewidth':1, 'linestyle':'--'})

#    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    d1 = scipy.stats.ranksums(data[0] , data[1])[1]
    d2 = scipy.stats.ranksums(data[0] , data[2])[1]
    d3 = scipy.stats.ranksums(data[1] , data[2])[1]

    
    ax.set_xticks([1 , 2 , 3])
    ax.set_xticklabels([label[0] , label[1] , label[2]] , fontsize = 28)
    ax.set_xlabel(label[0] + '_' + label[1] + ':' + str(d1) + ',' + label[0] + '_' + label[2] + ':' + str(d2) + ',' + label[1] + '_' + label[2] + ':'  + str(d3))
    ax.set_ylabel('log2(FPKM + 1)')
    ax.set_xlim((0.5 , 3.5))
    ax.set_ylim((vmin , vmax))
    return fig
    
def BA_AB_Box_plot(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'seagreen','linewidth':1},
            medianprops={'color':'seagreen','linewidth':1},
            capprops={'color':'seagreen','linewidth':1},
            whiskerprops={'color':'seagreen','linewidth':1, 'linestyle':'--'})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'chocolate','linewidth':1},
            medianprops={'color':'chocolate','linewidth':1},
            capprops={'color':'chocolate','linewidth':1},
            whiskerprops={'color':'chocolate','linewidth':1, 'linestyle':'--'})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7, 
            boxprops={'color': 'slateblue','linewidth':1},
            medianprops={'color':'slateblue','linewidth':1},
            capprops={'color':'slateblue','linewidth':1},
            whiskerprops={'color':'slateblue','linewidth':1, 'linestyle':'--'})

#    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    d1 = scipy.stats.ranksums(data[0] , data[1])[1]
    d2 = scipy.stats.ranksums(data[0] , data[2])[1]
    d3 = scipy.stats.ranksums(data[1] , data[2])[1]
    
    ax.set_xticks([1 , 2 , 3])
    ax.set_xticklabels(['IPSC_B_A' , 'IPSC_A_B' , 'stable'] , fontsize = 28)
    ax.set_xlabel('BA_AB:' + str(d1) + ',BA_stable:' + str(d2) + ',AB_stable:' + str(d3))
    ax.set_xlim((0.5 , 3.5))
#    ax.set_ylim((0 , 40))
    return fig
        
    
def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
  
    
pc_type = np.dtype({'names':['chr' , 'pc'] , 
                    'formats':['U8' , np.float]})
loc_type = np.dtype({'names':['chr' , 'pos'] , 
                    'formats':['U8' , np.int]})
lenth_type = np.dtype({'names':['chr' , 'len'] , 
                        'formats':['U8' , np.int]})
pc_type_1 = np.dtype({'names':['chr' , 'pos' , 'pc'] , 
                        'formats':['U8' , np.int , np.float]})

data_type_1 = np.dtype({'names':['chr' , 'pos' , 'CCs' , 'NTs' , 'fESC'] , 
                        'formats':['U8' , np.int , np.float , np.float , np.float]})
data_type_2 = np.dtype({'names':['chr' , 'pos' , 'MEF' , 'IPSC' , 'E14'] , 
                        'formats':['U8' , np.int , np.float , np.float , np.float]})


          
union_type_1 = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'F35' , 'F40' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
union_type_2 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'F35' , 'F40' , 'MEF' , 'IPS_P3' , 'E14'],
                 'formats':['U64' , 'U8' , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
  


chroms = ['1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , '18' , '19']
res = 200000


    
##Files

CCS = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\CCS_Traditonal_PC_200K_Compartment_200K.txt')
NT5 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NT5_Traditonal_PC_200K_Compartment_200K.txt')
NT6 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NT6_Traditonal_PC_200K_Compartment_200K.txt')
F35 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\F35_Traditonal_PC_200K_Compartment_200K.txt')
F40 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\F40_Traditonal_PC_200K_Compartment_200K.txt')
NTs = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\NTs_Traditonal_PC_200K_Compartment_200K.txt')
fESC = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Compartment_cooler\\compartment_new\\fESC_Traditonal_PC_200K_Compartment_200K.txt')




MEF = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\MEF_compartment_200K.txt')
IPS_P3 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P3_compartment_200K.txt')
IPS_P20 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPS_P20_compartment_200K.txt')
E14 = get_PC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\E14_compartment_200K.txt')






fESC_E14_common = []

for g in chroms:
    tmp_fESC = fESC[fESC['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    
    for i in range(len(tmp_fESC)):
        if tmp_fESC[i]['pos'] != tmp_E14[i]['pos']:
            print (g , i)
        if (tmp_fESC[i]['pc'] * tmp_E14[i]['pc'] > 0 ):
            fESC_E14_common.append((g , i * 200000 , (tmp_fESC[i]['pc'] + tmp_E14[i]['pc']) / 2))
fESC_E14_common = np.array(fESC_E14_common , dtype = pc_type_1)



    

#Nuclear transfer

##NT_process
A_B = []
B_A = []


for g in chroms:
    tmp_CCS = CCS[CCS['chr'] == g]
    tmp_NTs = NTs[NTs['chr'] == g]
    tmp_fESC = fESC[fESC['chr'] == g]
    tmp_common = fESC_E14_common[fESC_E14_common['chr'] == g]
    for i in range(len(tmp_CCS)):
        if (tmp_CCS[i]['pc'] > 0) and (tmp_fESC[i]['pc'] < 0):
            if i * 200000 in tmp_common['pos']:
                A_B.append((g , i * 200000 , tmp_CCS[i]['pc'] , tmp_NTs[i]['pc'] , tmp_fESC[i]['pc']))
        elif (tmp_CCS[i]['pc'] < 0) and (tmp_fESC[i]['pc'] > 0):
            if i * 200000 in tmp_common['pos']:
                B_A.append((g , i * 200000 , tmp_CCS[i]['pc'] , tmp_NTs[i]['pc'] , tmp_fESC[i]['pc']))
A_B = np.array(A_B , dtype = data_type_1)
B_A = np.array(B_A , dtype = data_type_1)




##IPSC

IPSC_A_B = []
IPSC_B_A = []
IPSC_stable = []
for g in chroms:
    tmp_MEF = MEF[MEF['chr'] == g]
    tmp_IPS_P3 = IPS_P3[IPS_P3['chr'] == g]
    tmp_E14 = E14[E14['chr'] == g]
    tmp_common = fESC_E14_common[fESC_E14_common['chr'] == g]
    for i in range(len(tmp_MEF)):
        if (tmp_MEF[i]['pc'] > 0) and (tmp_E14[i]['pc'] < 0):
            if i * 200000 in tmp_common['pos']:
                IPSC_A_B.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))
        elif (tmp_MEF[i]['pc'] < 0) and (tmp_E14[i]['pc'] > 0):
            if i * 200000 in tmp_common['pos']:
                IPSC_B_A.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))
        elif((tmp_MEF[i]['pc'] > 0) and (tmp_E14[i]['pc'] > 0)) or ((tmp_MEF[i]['pc'] < 0) and (tmp_E14[i]['pc'] < 0)):
            IPSC_stable.append((g , i * 200000 , tmp_MEF[i]['pc'] , tmp_IPS_P3[i]['pc'] , tmp_E14[i]['pc']))
IPSC_A_B = np.array(IPSC_A_B , dtype = data_type_2)
IPSC_B_A = np.array(IPSC_B_A , dtype = data_type_2)
IPSC_stable = np.array(IPSC_stable , dtype = data_type_2)
    
       

mm10 = np.loadtxt('E:\\Data\\literature_data\\genome\\mm10_chr.txt' , dtype = lenth_type)
gtf = Load_gtf('E:\\Data\\literature_data\\genome\\gencode.vM15.chr_patch_hapl_scaff.annotation.gtf')
# union_gene = get_union_gene_all_site_IPSC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC_ATAC\\RNA\\gene_expression\\NT_IPSC_all_genes.txt')
union_gene = get_union_gene_all_site_IPSC('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\IPSC\\RNA\\NT_IPSC_all_gene_expression.txt')



NT_c1 = [] ; NT_c2 = []
for i in A_B:
    chro = i[0]
    start = i[1]
    end = i[1] + 200000
    tmp = union_gene[union_gene['chr'] == 'chr' + chro]
    mask = (tmp['pos'] >= start) & (tmp['pos'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            NT_c1.append(j)
NT_c1 = np.array(NT_c1 , dtype = union_gene.dtype)

  

for i in B_A:
    chro = i[0]
    start = i[1]
    end = i[1] + 200000
    tmp = union_gene[union_gene['chr'] == 'chr' + chro]
    mask = (tmp['pos'] >= start) & (tmp['pos'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            NT_c2.append(j)
NT_c2 = np.array(NT_c2 , dtype = union_gene.dtype)                 



c1 = [] ; c2 = [] ; c3 = []
for i in IPSC_A_B:
    chro = i[0]
    start = i[1]
    end = i[1] + 200000
    tmp = union_gene[union_gene['chr'] == 'chr' + chro]
    mask = (tmp['pos'] >= start) & (tmp['pos'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c1.append(j)
c1 = np.array(c1 , dtype = union_gene.dtype)

            
for i in IPSC_B_A:
    chro = i[0]
    start = i[1]
    end = i[1] + 200000
    tmp = union_gene[union_gene['chr'] == 'chr' + chro]
    mask = (tmp['pos'] >= start) & (tmp['pos'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c2.append(j)
c2 = np.array(c2 , dtype = union_gene.dtype)


for i in IPSC_stable:
    chro = i[0]
    start = i[1]
    end = i[1] + 200000
    tmp = union_gene[union_gene['chr'] == 'chr' + chro]
    mask = (tmp['pos'] >= start) & (tmp['pos'] <= end)
    overlap = tmp[mask]
    if overlap.size > 0:
        for j in overlap:
            c3.append(j)
c3 = np.array(c3 , dtype = union_gene.dtype)




fig1 = Box_plot([np.log2(NT_c1['CCS'] + 1) , np.log2((NT_c1['NT5'] + NT_c1['NT6']) / 2 + 1) , np.log2((NT_c1['F35'] + NT_c1['F40']) / 2 + 1)] , ['CCS' , 'NTs' , 'fESC'] , -0.1 , 8)
fig2 = Box_plot([np.log2(NT_c2['CCS'] + 1) , np.log2((NT_c2['NT5'] + NT_c2['NT6']) / 2 + 1) , np.log2((NT_c2['F35'] + NT_c2['F40']) / 2 + 1)] , ['CCS' , 'NTs' , 'fESC'] , -0.1 , 8)

fig3 = Box_plot([np.log2(c1['MEF'] + 1) , np.log2(c1['IPSC'] + 1) , np.log2(c1['E14'] + 1)] , ['MEF' , 'iPSC' , 'E14'] , 2.5 , 4.5)
fig4 = Box_plot([np.log2(c2['MEF'] + 1) , np.log2(c2['IPSC'] + 1) , np.log2(c2['E14'] + 1)] , ['MEF' , 'iPSC' , 'E14'] , 2.5 , 4.5)





run_Plot(fig1 , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S5\\B_NT_ESC_compartment_classify_gene_FPKM_boxplot_A_B_2.pdf')
run_Plot(fig2 , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S5\\B_NT_ESC_compartment_classify_gene_FPKM_boxplot_B_A_2.pdf')
run_Plot(fig3 , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S5\\B_IPSC_ESC_compartment_classify_gene_FPKM_boxplot_A_B_2.pdf')
run_Plot(fig4 , 'F:\\work\\ntESC_3Dreprogramming\\Figures\\Plot_new\\S5\\B_IPSC_ESC_compartment_classify_gene_FPKM_boxplot_B_A_2.pdf')

out1 = open('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_compartment_A_B_geneNames.txt' , 'w')
out2 = open('H:\\Workspace_New\\data\\IPSC\\HiC\\Compartment\\compartment_new\\IPSC_compartment_A_B_geneID.txt' , 'w')

for i in c1:
    gene_name = i['gene_name']
    gene_id = gtf[gtf['gene_name'] == gene_name][0]['gene_id']
    out1.writelines(gene_name + '\n')
    out2.writelines(gene_id + '\n')
out1.close()
out2.close()



a = c1['IPSC'] - c1['MEF']
b = c2['IPSC'] - c2['MEF']
c = c3['IPSC'] - c3['MEF']

BA_AB_Box_plot([b,a,c])

