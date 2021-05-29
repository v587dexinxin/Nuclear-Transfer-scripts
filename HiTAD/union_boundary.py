# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 09:52:30 2020

@author: xxli
"""


from __future__ import division
import numpy as np
import os
import scipy
from scipy import stats
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles


tad_type = np.dtype({'names':['chr' , 'start' , 'end'],
                     'formats':['U4' , np.int , np.int]})

boun_type = np.dtype({'names':['chr' , 'site' ],
                     'formats':['U4' , np.int ]})

union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float]})
         
                 
                 

def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()
    
def Load_Tads(TadFil):
    Tads = np.loadtxt(TadFil , dtype = tad_type)
    return Tads


def Tads2Boundary(Tads):
    boundary =[]
    for i in Tads:
        chro = i['chr']
        start = i['start']
        end = i['end']
        boundary.append((chro , start))
        boundary.append((chro , end))
    boundary = list(set(boundary))
    boundary.sort(key = lambda x:(x[0] , x[1]))
    boundary = np.array(boundary , dtype = boun_type)
    return boundary
    
    
        
def union_boundary(data_list):
    boundary_all = []
    union = []
    for c in data_list:
        for i in c:
            boundary_all.append(i)          
    boundary_all.sort(key = lambda x:(x[0] , x[1]))
    
    boun_number = len(boundary_all)
    for i in range(boun_number-1 , -1 , -1):
        site = boundary_all[i]['site']
        start = boundary_all[i-1]['site'] - 40000
        end = boundary_all[i-1]['site'] + 40000
        if (site >= start) and (site <= end):
            boundary_all.remove(boundary_all[i])
        else:
            pass

    union = np.array(boundary_all , dtype =  boun_type)
    

    return union
    
    
    
def Common_boundary(boun_data , union_boundary):
    commom = []    
    
    for i in union_boundary:
        n = 0
        chro = i['chr']
        start = i['site'] - 40000
        end = i['site'] + 40000
        for j in boun_data:
            tmp = j[j['chr'] == chro]
            mask = (tmp['site'] >= start) & (tmp['site'] <= end)
            overlap = tmp[mask]
            if overlap.size != 0:
                n += 1
            else:
                pass
        if n == len(boun_data):
            commom.append(i)
    commom = np.array(commom , dtype = boun_type)
    return commom

def Common2_boundary(boun1 , boun2):
    common = []
    remain1 = []
    remain2 = []
    n1 = 0
    n2 = 0
#    len1 = len(Boun1)
#    len2 = len(Boun2)
#    if len1 < len2:
#        boun1 = Boun1
#        boun2 = Boun2
#    else:
#        boun1 = Boun2
#        boun2 = Boun1
        
    for i in boun1:
        chro = i['chr']
        start = i['site'] - 40000
        end = i['site'] + 40000
        tmp = boun2[boun2['chr'] == chro]
        mask = (tmp['site'] >= start) & (tmp['site'] <= end)
        overlap = tmp[mask]
        if overlap.size != 0:
            common.append(i)
            n1 += 1
        else:
            remain1.append(i)
    for i in boun2:
        chro = i['chr']
        start = i['site'] - 40000
        end = i['site'] + 40000
        tmp = boun1[boun1['chr'] == chro]
        mask = (tmp['site'] >= start) & (tmp['site'] <= end)
        overlap = tmp[mask]
        if overlap.size != 0:
            n2 += 1
        else:
            remain2.append(i)
    common = np.array(common , dtype = boun_type)
    remain1 = np.array(remain1 , dtype = boun_type)
    remain2 = np.array(remain2 , dtype = boun_type)
    
    return common , remain1 , remain2


def diff_boundary(boun1 , boun2):
    diff1 = []
    common = []
    for i in boun1:
        chro = i['chr']
        start = i['site'] - 40000
        end = i['site'] + 40000

        tmp = boun2[boun2['chr'] == chro]
        mask = (tmp['site'] >= start) & (tmp['site'] <= end) 
        overlap = tmp[mask]
        if overlap.size == 0:
            diff1.append(i)
        else:
            common.append(i)
    diff1 = np.array(diff1 , dtype = boun_type)
    common = np.array(common , dtype = boun_type)
    return diff1 , common

    
def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
 
def Box_plot(data):                
    left, bottom, width, height = 0.2 , 0.2 , 0.6 , 0.7
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.boxplot(data[0] , positions=[1] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'seagreen','linewidth':2},
            medianprops={'color':'seagreen','linewidth':2},
            capprops={'color':'seagreen','linewidth':2},
            whiskerprops={'color':'seagreen','linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[1] , positions=[2] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'chocolate','linewidth':2},
            medianprops={'color':'chocolate','linewidth':2},
            capprops={'color':'chocolate','linewidth':2},
            whiskerprops={'color':'chocolate','linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[2] , positions=[3] , showfliers=False, widths = 0.7, 
            boxprops={'color': 'slateblue','linewidth':2},
            medianprops={'color':'slateblue','linewidth':2},
            capprops={'color':'slateblue','linewidth':2},
            whiskerprops={'color':'slateblue','linewidth':2, 'linestyle':'--'})
    ax.boxplot(data[3] , positions=[4] , showfliers=False, widths = 0.7 ,
            boxprops={'color': 'deeppink','linewidth':2},
            medianprops={'color':'deeppink','linewidth':2},
            capprops={'color':'deeppink','linewidth':2},
            whiskerprops={'color':'deeppink','linewidth':2, 'linestyle':'--'})
#    ax.plot([0.5,4.5],[0,0], lw = 1.5, ls = '--', color = 'darkblue')
    d1 = scipy.stats.ranksums(data[0] , data[1])[1]
    d2 = scipy.stats.ranksums(data[0] , data[2])[1]
    d3 = scipy.stats.ranksums(data[0] , data[3])[1]
    d4 = scipy.stats.ranksums(data[1] , data[2])[1]
    d5 = scipy.stats.ranksums(data[1] , data[3])[1]
    d6 = scipy.stats.ranksums(data[2] , data[3])[1]
    
    ax.set_xticks([1 , 2 , 3 , 4])
    ax.set_xticklabels(['CCs' , 'NT5' , 'NT6' , 'F40' ] , fontsize = 28)
    ax.set_xlabel('CCS_NT5:' + str(d1) + ',CCS_NT6:' + str(d2) + ',CCS_fESC:' + str(d3) + '\n' + 'NT5_NT6:' + str(d4) + ',NT5_fESC:' + str(d5) + ',NT6_fESC:' + str(d6))
    ax.set_xlim((0.5 , 4.5))
    ax.set_ylim((0 , 20))
    return fig

def Sig_To_1K(fil):
    """
    """
    sig_type = np.dtype({'names':['chr','start' , 'end' , 'score'],
                      'formats':['S4',np.int , np.int , np.float]})
    RNAData = np.loadtxt(fil , dtype=sig_type)
    
    chroms = set(RNAData['chr'])
    New_Data = {}
    for c in chroms:
        New_Data[c] = {}
        tmp_data = RNAData[RNAData['chr'] == c]
        max_ = tmp_data['end'].max()
        bin_size = max_ // 1000 + 1
        New_Data[c] = np.zeros((bin_size,))
        for line in tmp_data:
            start = line['start'] // 1000
#            end = line['end'] // 1000
#            for i in range(start,end):
#                New_Data[c][i] += line['score'] /(end - start)
            New_Data[c][start] += line['score'] 
    
    return New_Data
       
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    
def Write2genes(filname , gene):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['Gene_name' , 'Chr' , 'Gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')
        for i in gene:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()  


def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()

                      
            
TadFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain'
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19']
#cells = ['NT5','NT6','fESC']

CCS = Load_Tads(os.path.join(TadFolder , 'CCS_Domain_bottom_40K_respective_stable_400K.txt'))
NT5 = Load_Tads(os.path.join(TadFolder , 'NT5_Domain_bottom_40K_respective_stable_400K.txt'))
NT6 = Load_Tads(os.path.join(TadFolder , 'NT6_Domain_bottom_40K_respective_stable_400K.txt'))
F35 = Load_Tads(os.path.join(TadFolder , 'F35_Domain_bottom_40K_respective_stable_400K.txt'))
F40 = Load_Tads(os.path.join(TadFolder , 'F40_Domain_bottom_40K_respective_stable_400K.txt'))




CCS_b_0 = Tads2Boundary(CCS)
NT5_b_0 = Tads2Boundary(NT5)
NT6_b_0 = Tads2Boundary(NT6)
F35_b_0 = Tads2Boundary(F35)
F40_b_0 = Tads2Boundary(F40)


CCS_b = union_boundary([CCS_b_0])
NT5_b = union_boundary([NT5_b_0])
NT6_b = union_boundary([NT6_b_0])
F35_b = union_boundary([F35_b_0])
F40_b = union_boundary([F40_b_0])


common_NT_b , n1 , n2= Common2_boundary(NT5_b , NT6_b)
common_fESC_b , n1 , n2 = Common2_boundary(F35_b , F40_b)
union_CCS_NT_fESC = union_boundary([CCS_b , common_NT_b , commom_fESC_b])

##union_boundaries
Write2fils_nochr(os.path.join(TadFolder , 'union_Filtered_Boundary_40K.txt') , union_CCS_NT_fESC) 


    