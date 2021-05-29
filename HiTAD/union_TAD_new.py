# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 10:20:24 2021

@author: xxli
"""

from __future__ import division
import numpy as np
#from tadlib.calfea.analyze import getmatrix
import matplotlib
# Use a non-interactive backend
# matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import os
from scipy.special import ndtr
import math
import csv
from itertools import islice
import pandas as pd

#--------------------------------------------------------------------------
## Matplotlib Settings
matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib_venn import venn2, venn2_circles
from matplotlib_venn import venn3, venn3_circles
# Our Own Color Map
my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
my_cmap.set_bad('#2672a1')




def Get_tads(LoopFil):
    """
    cell : String ,The name of cells  e.g.'fESC','ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['U8' , np.int , np.int]})
    loops = []
    fileSource = os.path.join(LoopFil)
    Loop = np.loadtxt(fileSource, usecols = (0,1,2) , dtype = loop_type)
    Loop = Loop[Loop['chr'] != 'X']
    
    for i in Loop:
        if i['end'] - i['start'] >= 400000:
            loops.append((i['chr'].lstrip('chr') , i['start'] , i['end']))
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops


def Get_union_tads(Loop_list , cell_line , initial_distance , raw = True):
    loop_type_1 = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['U8' , np.int , np.int , 'U8']})
    all_loops = []
    for i in range(len(cell_line)):
        for j in Loop_list[i]:
            all_loops.append((j['chr'] , j['start'] , j['end'] , cell_line[i]))
            
    all_loops = np.array(all_loops , dtype = loop_type_1)
    union = clustering(all_loops , initial_distance)
    union_new = clustering_rescue(union , R , raw)
    return union_new

        


def Sort(a , s_1 , s_2):
    '''
    a: list of needed to be sorted
    '''
    a.sort(key = lambda x:(x[s_1],x[s_2]))
    return a



def center_sites(lists):
    sum_x = 0
    sum_y = 0
    for x in lists:
        sum_x += x[1]
        sum_y += x[2]
    n = len(lists)
    return [float(sum_x)/n, float(sum_y)/n]

def distance(sites_1,sites_2):
    dis = math.sqrt((sites_1[1] / 40000 - sites_2[1] / 40000)**2 + (sites_1[2] / 40000-sites_2[2] / 40000)**2)
    return dis

    
def clustering(loop_sites, dis):
    chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19']
    classes = []
    for i in chrom:
        c_loops = sorted(loop_sites[loop_sites['chr'] == i], key = lambda x:x[1])
        while True :
            c_classs = []
            c_classs.append((c_loops[0][0] , c_loops[0][1] , c_loops[0][2] , c_loops[0][3]))
            c_loops.remove(c_loops[0])
            center = center_sites(c_classs)
            for loop in c_loops:     
                 if distance([i , center[0] , center[1]], loop) <= dis:
                      c_classs.append((loop[0] , loop[1] , loop[2] , loop[3]))
                      center = center_sites(c_classs)
                      c_loops.remove(loop)
            classes.append(c_classs)
            if len(c_loops) == 0:
                break
    return classes

def clustering_rescue(clustered_sites , dis , raw = True):
    loop_type_1 = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['U8' , np.int , np.int , 'U8']})
    single = []
    multi = []
    for i in clustered_sites:
        if len(i) > 1:
            multi.append([s for s in i])
        else:
            single.append((i[0][0] , i[0][1] , i[0][2] , i[0][3]))
        
    single = np.array(single , dtype = loop_type_1)
    for i in range(len(multi)):
        for j in multi[i]:
            chro = j[0]
            s_s = j[1] - dis
            s_e = j[1] + dis
            e_s = j[2] - dis
            e_e = j[2] + dis
            s = single[single['chr'] == chro]
            mask = (s['start'] >= s_s) & (s['start'] <= s_e) & (s['end'] >= e_s) & (s['end'] <= e_e)
            overlap = s[mask]
            if overlap.size != 0:
                single = list(single)
                for k in overlap:
                    clustered_sites[clustered_sites.index(multi[i])].append((k[0] , k[1] , k[2] , k[3]))
                    multi[i].append((k[0] , k[1] , k[2] , k[3]))
                    clustered_sites.remove([(k[0] , k[1] , k[2] , k[3])])
                    single.remove(k)
                single = np.array(single , dtype = loop_type_1)
            else:
                continue
                

    union_new = [np.array(a , dtype = loop_type_1) for a in clustered_sites]
    union_sites = []
    for i in union_new:
        chro = i[0][0]
        center = center_sites(i)
        union_sites.append((chro , int(center[0]) , int(center[1])))
    union_sites = Sort(union_sites , 0 , 1)
    union_sites = np.array(union_sites , dtype = tad_type)
    
    if raw == True:
        return union_new
    else:
        return union_sites


def Get_center_sites(site_list):
    sites = []
    for i in site_list:
        chro = i[0][0]
        center = center_sites(i)
        sites.append((chro , int(center[0]) , int(center[1])))
    sites = Sort(sites , 0 , 1)
    sites = np.array(sites , dtype = tad_type)
    return sites

def Get_diff_loops(loop_1 , loop_2):
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['U8' , np.int , np.int]})
    common_1 = []
    common_2 = []
    for i in loop_1:
        chro = i['chr']
        tmp = loop_2[loop_2['chr'] == chro]
        for j in tmp:
            d = distance([chro , i[1] , i[2]] , j)
            if d <= initial_distance:
                common_1.append(i)
    for i in loop_2:
        chro = i['chr']
        tmp = loop_1[loop_1['chr'] == chro]
        for j in tmp:
            d = distance([chro , i[1] , i[2]] , j)
            if d <= initial_distance:
                common_2.append(i)
    common_1 = np.array(common_1 , dtype = loop_type)
    common_2 = np.array(common_2 , dtype = loop_type)
    loop_1 = list(loop_1)
    loop_2 = list(loop_2)
    print (len(loop_1) , len(loop_2) , len(common_1) , len(common_2))
    for i in common_1:
        try:
            loop_1.remove(i)
        except:
            pass
    for j in common_2:
        try:
            loop_2.remove(j)
        except:
            pass
    print (len(loop_1) , len(loop_2))   
    diff_1 = np.array(Sort(loop_1 , 0 , 1 ), dtype = loop_type)
    diff_2 = np.array(Sort(loop_2 , 0 , 1 ), dtype = loop_type)
    return diff_1 , diff_2 , common_1 , common_2


def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
        
        
        
def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            out.writelines('\t'.join([str(x) for x in i]) +'\n')
    out.close()



def Write2genes(filname , gene):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['Gene_name' , 'Chr' , 'Gene_site' , 'CCS' , 'NT5' , 'NT6' , 'F35' , 'F40']) + '\n')
        for i in gene:
            out.writelines('\t'.join([i['gene_name'] , i['chr'] , str(i['gene_site']) , str(i['CCS_FPKM']) , str(i['NT5_FPKM']) , str(i['NT6_FPKM']) , str(i['F35_FPKM']) , str(i['F40_FPKM']) ])+'\n')
    out.close()          
            
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()



def get_union_gene_sites(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8' , 'U8' , np.int , np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
    union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})
         
                 
    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    union = []
    for i in union_gene:
        # if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
            gene_site = int((i['start'] + i['end']) / 2)
            union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS_FPKM'] , i['NT2_FPKM'] , i['NT3_FPKM'] , i['NT4_FPKM'] , i['NT5_FPKM'] , i['NT6_FPKM'] , i['F35_FPKM'] , i['F37_FPKM'] , i['F40_FPKM'] , i['F41_FPKM']))
    union = np.array(union , dtype = union_type_1)
    return union


def Load_gtf(gtfil):
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



tad_type = ({'names':['chr' , 'start' , 'end'],
              'formats':['U8' , np.int , np.int]})
              
union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float]})

union_type_1 = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS_FPKM' , 'NT2_FPKM' , 'NT3_FPKM' , 'NT4_FPKM' , 'NT5_FPKM' , 'NT6_FPKM' , 'F35_FPKM' , 'F37_FPKM' , 'F40_FPKM' , 'F41_FPKM'],
                 'formats':['U64' , 'U8', np.int , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float , np.float]})

         
                 
R = 40000
initial_distance = math.sqrt(2) + 0.05


TadFolder = 'F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\'        
CCSFil = 'CCS_Domain_bottom_40K_respective_stable_400K.txt'
NT5Fil = 'NT5_Domain_bottom_40K_respective_stable_400K.txt'
NT6Fil = 'NT6_Domain_bottom_40K_respective_stable_400K.txt'
F35Fil = 'F35_Domain_bottom_40K_respective_stable_400K.txt'
F40Fil = 'F40_Domain_bottom_40K_respective_stable_400K.txt'
NTsFil = 'NTs_Domain_bottom_40K_respective_stable_400K.txt'
fESCFil = 'fESC_Domain_bottom_40K_respective_stable_400K.txt'


CCS_tad = Get_tads(os.path.join(TadFolder , CCSFil))
NT5_tad = Get_tads(os.path.join(TadFolder , NT5Fil))
NT6_tad = Get_tads(os.path.join(TadFolder , NT6Fil))
F35_tad = Get_tads(os.path.join(TadFolder , F35Fil))
F40_tad = Get_tads(os.path.join(TadFolder , F40Fil))
NTs_tad = Get_tads(os.path.join(TadFolder , NTsFil))
fESC_tad = Get_tads(os.path.join(TadFolder , fESCFil))



##-----------------------------respective_common---------------------------------------------------
NT5_1 , NT6_1 , common_NT , common_NT_2 = Get_diff_loops(NT5_tad , NT6_tad)
F35_1 , F40_1 , common_fESC , common_fESC_2 = Get_diff_loops(F35_tad , F40_tad)


union = Get_union_tads([CCS_tad , common_NT , common_fESC] , ['CCS' , 'NT' , 'fESC'] , initial_distance , True)

union_tads = Get_center_sites(union)


Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\union_domain_NTs_fESC_common.txt' , union_tads)



##-----------------------------merged---------------------------------------------------



union = Get_union_tads([CCS_tad , NTs_tad , fESC_tad] , ['CCS' , 'NT' , 'fESC'] , initial_distance , True)

union_tads = Get_center_sites(union)


Write2fils_nochr('F:\\work\\ntESC_3Dreprogramming\\Workspace_New\\data\\HiC\\Domain_new\\xtwang\\bottom_domain\\union_domain_NTs_fESC_merged.txt' , union_tads)






















