# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 17:04:50 2020

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
import csv

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




def Get_loops(LoopFil):
    """
    cell : String ,The name of cells  e.g.'fESC','ESC' , 'NT5' , 'NT6' , 'CCS'
    get the loops' location as a dictionary
    """
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['S8' , np.int , np.int]})
    loops = []
    fileSource = os.path.join(LoopFil)
    Loop = np.loadtxt(fileSource, usecols = (0,1,2) , dtype = loop_type, skiprows = 1)
    Loop = Loop[Loop['chr'] == 'X']
    
    for i in Loop:
        if i['end'] - i['start'] >= 300000:
            loops.append(i)
        else:
            continue
    loops = np.array(loops , dtype = loop_type)
    return loops

def Get_union_loops(Loop_list , cell_line , initial_distance , raw = True):
    loop_type_1 = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['S8' , np.int , np.int , 'S8']})
    all_loops = []
    for i in range(len(cell_line)):
        for j in Loop_list[i]:
            all_loops.append((j['chr'] , j['start'] , j['end'] , cell_line[i]))
            
    all_loops = np.array(all_loops , dtype = loop_type_1)
    union = clustering(all_loops , initial_distance)
    union_new = clustering_rescue(union , 2 * R , raw)
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
    dis = math.sqrt((sites_1[0] / 20000 - sites_2[1] / 20000)**2 + (sites_1[1] / 20000-sites_2[2] / 20000)**2)
    return dis

    
def clustering(loop_sites, dis):
    chrom = ['X']
    classes = []
    for i in chrom:
        c_loops = sorted(loop_sites[loop_sites['chr'] == i], key = lambda x:x[1])
        while True :
            c_classs = []
            c_classs.append((c_loops[0][0] , c_loops[0][1] , c_loops[0][2] , c_loops[0][3]))
            c_loops.remove(c_loops[0])
            center = center_sites(c_classs)
            for loop in c_loops:     
                 if distance(center, loop) <= dis:
                      c_classs.append((loop[0] , loop[1] , loop[2] , loop[3]))
                      center = center_sites(c_classs)
                      c_loops.remove(loop)
            classes.append(c_classs)
            if len(c_loops) == 0:
                break
    return classes

def clustering_rescue(clustered_sites , dis , raw = True):
    loop_type_1 = ({'names':['chr' , 'start' , 'end' , 'cell'],
                  'formats':['S8' , np.int , np.int , 'S8']})
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
    union_sites = np.array(union_sites , dtype = loop_type)
    
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
    sites = np.array(sites , dtype = loop_type)
    return sites

def Get_diff_loops(loop_1 , loop_2):
    loop_type = ({'names':['chr' , 'start' , 'end'],
                  'formats':['S8' , np.int , np.int]})
    common_1 = []
    common_2 = []
    for i in loop_1:
        chro = i['chr']
        tmp = loop_2[loop_2['chr'] == chro]
        for j in tmp:
            d = distance([i[1] , i[2]] , j)
            if d <= initial_distance:
                common_1.append(i)
    for i in loop_2:
        chro = i['chr']
        tmp = loop_1[loop_1['chr'] == chro]
        for j in tmp:
            d = distance([i[1] , i[2]] , j)
            if d <= initial_distance:
                common_2.append(i)
    common_1 = np.array(common_1 , dtype = loop_type)
    common_2 = np.array(common_2 , dtype = loop_type)
    loop_1 = list(loop_1)
    loop_2 = list(loop_2)
    print len(loop_1) ; print len(loop_2) ; print len(common_1) ; print len(common_2)
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
    print len(loop_1) ; print len(loop_2)   
    diff_1 = np.array(Sort(loop_1 , 0 , 1 ), dtype = loop_type)
    diff_2 = np.array(Sort(loop_2 , 0 , 1 ), dtype = loop_type)
    return diff_1 , diff_2 , common_1 , common_2


def get_union_gene(Fil):
    union_type = ({'names':['gene_name' , 'chr' , 'strand' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8' , 'S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

    union_gene = np.loadtxt(Fil , dtype = union_type , skiprows = 1)
    return union_gene
        
        
        
def Write2fils_nochr(filname , peaks):
    with open(filname,'w') as out:
        for i in peaks:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()            

def Write2genes(filname , gene):
    with open(filname,'w') as out:
        out.writelines('\t'.join(['Gene_name' , 'Chr' , 'Gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC']) + '\n')
        for i in gene:
            i = np.array(list(i),dtype = str)
            out.writelines('\t'.join(i)+'\n')
    out.close()                        
            
def Get_common_loops(union):
    common = []
    for i in union:
        if len(i) == 3:
            common.append(i)
    common = Get_center_sites(common)
    return common
            
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
            

loop_type = ({'names':['chr' , 'start' , 'end'],
              'formats':['S8' , np.int , np.int]})
              
union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})
         
                 
R = 20000
initial_distance = 2 * math.sqrt(2) + 0.05


LoopFolder = 'H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05'       
CCSFil = 'Cluster_CCS_loops20K_3.txt'
NT5Fil = 'Cluster_NT5_loops20K_3.txt'
NT6Fil = 'Cluster_NT6_loops20K_3.txt'
fESCFil = 'Cluster_fESC_loops20K_3.txt'



CCS_loop = Get_loops(os.path.join(LoopFolder , CCSFil))
NT5_loop = Get_loops(os.path.join(LoopFolder , NT5Fil))
NT6_loop = Get_loops(os.path.join(LoopFolder , NT6Fil))
fESC_loop = Get_loops(os.path.join(LoopFolder , fESCFil))

###------------------------------Loop_Venn4------------------------------------------


CCS_1 , NT5_1 , common_12 , common_12_2 = Get_diff_loops(CCS_loop , NT5_loop)
CCS_1 , NT6_1 , common_13 , common_13_2 = Get_diff_loops(CCS_loop , NT6_loop)
CCS_1 , fESC_1 , common_14 , common_14_2 = Get_diff_loops(CCS_loop , fESC_loop)
NT5_1 , NT6_1 , common_23 , common_23_2 = Get_diff_loops(NT5_loop , NT6_loop)
NT5_1 , fESC_1 , common_24 , common_24_2 = Get_diff_loops(NT5_loop , fESC_loop)
NT6_1 , fESC_1 , common_34 , common_34_2 = Get_diff_loops(NT6_loop , fESC_loop)

union_123 = Get_union_loops([CCS_loop , NT5_loop , NT6_loop] , ['CCS' , 'NT5' , 'NT6'] , initial_distance , True)
union_124 = Get_union_loops([CCS_loop , NT5_loop , fESC_loop] , ['CCS' , 'NT5' , 'fESC'] , initial_distance , True)
union_134 = Get_union_loops([CCS_loop , NT6_loop , fESC_loop] , ['CCS' , 'NT6' , 'fESC'] , initial_distance , True)
union_234 = Get_union_loops([NT5_loop , NT6_loop , fESC_loop] , ['NT5' , 'NT6' , 'fESC'] , initial_distance , True)
union_1234 = Get_union_loops([CCS_loop , NT5_loop , NT6_loop , fESC_loop] , ['CCS' , 'NT5' , 'NT6' , 'fESC'] , initial_distance , True)

common_123 = Get_common_loops(union_123)
common_124 = Get_common_loops(union_124)
common_134 = Get_common_loops(union_134)
common_234 = Get_common_loops(union_234)
common_1234 = Get_common_loops(union_1234)
        
print 'area1=' + str(len(CCS_loop))
print 'area2=' + str(len(NT5_loop))
print 'area3=' + str(len(NT6_loop))
print 'area4=' + str(len(fESC_loop))
print 'n12=' + str(len(common_12)) 
print 'n13=' + str(len(common_13))
print 'n14=' + str(len(common_14))
print 'n23=' + str(len(common_23))
print 'n24=' + str(len(common_24))
print 'n34=' + str(len(common_34))
print 'n123=' + str(len(common_123))
print 'n124=' + str(len(common_124))
print 'n134=' + str(len(common_134))
print 'n234=' + str(len(common_234))
print 'n1234=' + str(len(common_1234))



###----------------------------Loop_Venn3_CCS_NTs_fESC-------------------------------------------------
NT5_1 , NT6_1 , common_NT , common_NT_2 = Get_diff_loops(NT5_loop , NT6_loop)


union = Get_union_loops([CCS_loop , common_NT , fESC_loop] , ['CCS' , 'NT' , 'fESC'] , initial_distance , True)


Abc = 0 ; aBc = 0 ; ABc = 0 ; abC = 0 ; AbC = 0 ; aBC = 0 ; ABC = 0
c1  = [] ; c2 = [] ; c3 = [] ; c4 = [] ;  c5 = [] ; c6 = [] ; c7 = []

for i in union:
    ce = i['cell']
    if ('CCS' in ce) and ('NT' not in ce) and ('fESC' not in ce):
        Abc += 1
        c1.append(i)
    if ('CCS' not in ce) and ('NT' in ce) and ('fESC' not in ce):
        aBc += 1
        c2.append(i)
    if ('CCS' in ce) and ('NT' in ce) and ('fESC' not in ce):
        ABc += 1
        c3.append(i)
    if ('CCS' not in ce) and ('NT' not in ce) and ('fESC' in ce):
        abC += 1
        c4.append(i)
    if ('CCS' in ce) and ('NT' not in ce) and ('fESC' in ce):
        AbC += 1
        c5.append(i)
    if ('CCS' not in ce) and ('NT' in ce) and ('fESC' in ce):
        aBC += 1
        c6.append(i)
    if ('CCS' in ce) and ('NT' in ce) and ('fESC' in ce):
        ABC += 1
        c7.append(i) 
        
print Abc , aBc , ABc , abC , AbC , aBC , ABC

fig = plt.figure(figsize = (10, 10))
venn3(subsets=(Abc , aBc , ABc , abC , AbC , aBC , ABC), set_labels=('CCS', 'NT' , 'fESC'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S8_figs\\CCS_NT(common)fESC_Loop_ChrX_Venn3.pdf')

loops = {'CCS_noNT_nofESC':Get_center_sites(c1) , 'noCCS_NT_nofESC':Get_center_sites(c2) , 'CCS_NT_nofESC' : Get_center_sites(c3),
         'noCCS_noNT_fESC':Get_center_sites(c4) , 'CCS_noNT_fESC':Get_center_sites(c5) , 'noCCS_NT_fESC' : Get_center_sites(c6) ,
         'CCS_NT_fESC' : Get_center_sites(c7)}

union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)


genes = {}

for cl in loops:
    genes[cl] = []
    for i in loops[cl]:
        chro = i['chr']
        start_s = i['start'] - 20000
        start_e = i['start'] + 20000
        end_s = i['end'] - 20000
        end_e = i['end'] + 20000
        tmp = union[union['chr'] == 'chr' + chro]
        mask1 = (tmp['gene_site'] >= start_s) & (tmp['gene_site'] <= start_e)
        mask2 = (tmp['gene_site'] >= end_s) & (tmp['gene_site'] <= end_e)
        overlap1 = tmp[mask1]
        overlap2 = tmp[mask2]
        if overlap1.size != 0:
            for j in overlap1:
                genes[cl].append(j)
        if overlap2.size != 0:
            for j in overlap2:
                genes[cl].append(j)
    genes[cl] = np.array(genes[cl] , dtype = union_type)
    
    
for k,v in genes.items():
    print k , len(v)
    Write2genes('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\' + k + '_genes_ChrX.txt' , v)                


##-----------------------------------------Loop_genes_Venn2-------------------------------------

loops_1 = np.hstack((loops['noCCS_NT_fESC'] , loops['CCS_NT_fESC']))     
genes = []

for i in loops_1:
    chro = i['chr']
    start_s = i['start'] - 20000
    start_e = i['start'] + 20000
    end_s = i['end'] - 20000
    end_e = i['end'] + 20000
    tmp = union[union['chr'] == 'chr' + chro]
    mask1 = (tmp['gene_site'] >= start_s) & (tmp['gene_site'] <= start_e)
    mask2 = (tmp['gene_site'] >= end_s) & (tmp['gene_site'] <= end_e)
    overlap1 = tmp[mask1]
    overlap2 = tmp[mask2]
    if overlap1.size != 0:
        for j in overlap1:
            genes.append(j)
    if overlap2.size != 0:
        for j in overlap2:
            genes.append(j)
            
genes = np.array(genes , dtype = union_type)

NT5_diff = csv.reader(open('H:\\Workspace_New\\data\\RNA\\diff_expression\\New_fc_1.5\\Filtered_NT5_fESC.csv' , 'r'))
NT6_diff = csv.reader(open('H:\\Workspace_New\\data\\RNA\\diff_expression\\New_fc_1.5\\Filtered_NT6_fESC.csv' , 'r'))

diff_gene = []
for i in NT5_diff:
    break
for i in NT6_diff:
    break
    
for i in [NT5_diff , NT6_diff]:
    for j in i:
        diff_gene.append(j[7])

diff_gene = set(diff_gene)      
  
diff = []        
for i in genes:
    if i['gene_name'] in diff_gene:
        diff .append(i)
diff = np.array(diff , dtype = union_type)

percent_l = len(diff) / len(genes)



###----------------------------Loop_Venn3_NT5_NT6_fESC-------------------------------------------------



union = Get_union_loops([NT5_loop , NT6_loop , fESC_loop] , ['NT5' , 'NT6' , 'fESC'] , initial_distance , True)


Abc = 0 ; aBc = 0 ; ABc = 0 ; abC = 0 ; AbC = 0 ; aBC = 0 ; ABC = 0
c1  = [] ; c2 = [] ; c3 = [] ; c4 = [] ;  c5 = [] ; c6 = [] ; c7 = []

for i in union:
    ce = i['cell']
    if ('NT5' in ce) and ('NT6' not in ce) and ('fESC' not in ce):
        Abc += 1
        c1.append(i)
    if ('NT5' not in ce) and ('NT6' in ce) and ('fESC' not in ce):
        aBc += 1
        c2.append(i)
    if ('NT5' in ce) and ('NT6' in ce) and ('fESC' not in ce):
        ABc += 1
        c3.append(i)
    if ('NT5' not in ce) and ('NT6' not in ce) and ('fESC' in ce):
        abC += 1
        c4.append(i)
    if ('NT5' in ce) and ('NT6' not in ce) and ('fESC' in ce):
        AbC += 1
        c5.append(i)
    if ('NT5' not in ce) and ('NT6' in ce) and ('fESC' in ce):
        aBC += 1
        c6.append(i)
    if ('NT5' in ce) and ('NT6' in ce) and ('fESC' in ce):
        ABC += 1
        c7.append(i) 
        
print Abc , aBc , ABc , abC , AbC , aBC , ABC

fig = plt.figure(figsize = (10, 10))
venn3(subsets=(Abc , aBc , ABc , abC , AbC , aBC , ABC), set_labels=('NT5', 'NT6' , 'fESC'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S8_figs\\NT5_NT6_fESC_Loop_ChrX_Venn3.pdf')

loops = {'NT5_noNT6_nofESC':Get_center_sites(c1) , 'noNT5_NT6_nofESC':Get_center_sites(c2) , 'NT5_NT6_nofESC' : Get_center_sites(c3),
         'noNT5_noNT6_fESC':Get_center_sites(c4) , 'NT5_noNT6_fESC':Get_center_sites(c5) , 'noNT5_NT6_fESC' : Get_center_sites(c6) ,
         'NT5_NT6_fESC' : Get_center_sites(c7)}

union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)


genes = {}

for cl in loops:
    genes[cl] = []
    for i in loops[cl]:
        chro = i['chr']
        start_s = i['start'] - 20000
        start_e = i['start'] + 20000
        end_s = i['end'] - 20000
        end_e = i['end'] + 20000
        tmp = union[union['chr'] == 'chr' + chro]
        mask1 = (tmp['gene_site'] >= start_s) & (tmp['gene_site'] <= start_e)
        mask2 = (tmp['gene_site'] >= end_s) & (tmp['gene_site'] <= end_e)
        overlap1 = tmp[mask1]
        overlap2 = tmp[mask2]
        if overlap1.size != 0:
            for j in overlap1:
                genes[cl].append(j)
        if overlap2.size != 0:
            for j in overlap2:
                genes[cl].append(j)
    genes[cl] = np.array(genes[cl] , dtype = union_type)
    
    
for k,v in genes.items():
    print k , len(v)
    Write2genes('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_ChrX_Venn3(NT5_NT6_fESC)\\' + k + '_genes_ChrX.txt' , v)                


