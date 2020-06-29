# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:50:04 2020

@author: han-luo
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
                     'formats':['S4' , np.int , np.int]})

boun_type = np.dtype({'names':['chr' , 'site' ],
                     'formats':['S4' , np.int ]})

union_type = ({'names':['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                 'formats':['S64' , 'S8', np.int , np.float , np.float , np.float , np.float]})
         
                 
                 

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
            
TadFolder = 'H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain'
chrom = ['1', '2', '3', '4', '5', '6', '7', '8', '9','10','11', '12', '13', '14', '15', '16', '17', '18', '19']
#cells = ['NT5','NT6','fESC']

CCS = Load_Tads(os.path.join(TadFolder , 'CCS_40K_allreps_union_small_domain.txt'))
CCS = CCS[CCS['chr'] == 'X']
NT5 = Load_Tads(os.path.join(TadFolder , 'NT5_40K_allreps_union_small_domain.txt'))
NT5 = NT5[NT5['chr'] == 'X']
NT6 = Load_Tads(os.path.join(TadFolder , 'NT6_40K_allreps_union_small_domain.txt'))
NT6 = NT6[NT6['chr'] == 'X']
fESC = Load_Tads(os.path.join(TadFolder , 'fESC_40K_allreps_union_small_domain.txt'))
fESC = fESC[fESC['chr'] == 'X']

CCS_b_0 = Tads2Boundary(CCS)
NT5_b_0 = Tads2Boundary(NT5)
NT6_b_0 = Tads2Boundary(NT6)
fESC_b_0 = Tads2Boundary(fESC)

CCS_b = union_boundary([CCS_b_0])
NT5_b = union_boundary([NT5_b_0])
NT6_b = union_boundary([NT6_b_0])
fESC_b = union_boundary([fESC_b_0])


###----------------------------Venn4--------------------------------
common_12 , n1 , n2= Common2_boundary(CCS_b , NT5_b)
common_13 , n1 , n2= Common2_boundary(CCS_b , NT6_b)
common_14 , n1 , n2= Common2_boundary(CCS_b , fESC_b)
common_23 , n1 , n2= Common2_boundary(NT5_b , NT6_b)
common_24 , n1 , n2= Common2_boundary(NT5_b , fESC_b)
common_34 , n1 , n2= Common2_boundary(NT6_b , fESC_b)


union_CCS_NT5_NT6 = union_boundary([CCS_b , NT5_b , NT6_b])
union_CCS_NT5_fESC = union_boundary([CCS_b , NT5_b , fESC_b])
union_CCS_NT6_fESC = union_boundary([CCS_b , NT6_b , fESC_b])
union_NT5_NT6_fESC = union_boundary([NT5_b , NT6_b , fESC_b])
union_all = union_boundary([CCS_b , NT5_b , NT6_b , fESC_b])



common_123 = Common_boundary([CCS_b , NT5_b , NT6_b] , union_CCS_NT5_NT6)
common_124 = Common_boundary([CCS_b , NT5_b , fESC_b] , union_CCS_NT5_fESC)
common_134 = Common_boundary([CCS_b , NT6_b , fESC_b] , union_CCS_NT6_fESC)
common_234 = Common_boundary([NT5_b , NT6_b , fESC_b] , union_NT5_NT6_fESC)
common_1234 = Common_boundary([CCS_b , NT5_b , NT6_b , fESC_b] , union_all)


print 'area1=' + str(len(CCS_b))
print 'area2=' + str(len(NT5_b))
print 'area3=' + str(len(NT6_b))
print 'area4=' + str(len(fESC_b))
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




###---------------------------Venn3_CCs_NTs_fESC----------------------------------------


common_NT_b , n1 , n2= Common2_boundary(NT5_b , NT6_b)
union_CCS_NT_fESC = union_boundary([CCS_b , common_NT_b , fESC_b])

CCS_NT , CCS_noNT , noCCS_NT = Common2_boundary(CCS_b , common_NT_b)
CCS_fESC , CCS_nofESC , noCCS_fESC  = Common2_boundary(CCS_b , fESC_b)
NT_fESC , NT_nofESC , noNT_fESC = Common2_boundary(common_NT_b , fESC_b)
CCS_NT_fESC = Common_boundary([CCS_b , common_NT_b , fESC_b] , union_CCS_NT_fESC)

CCS_noNT_nofESC , CCS_noNT_fESC= diff_boundary(CCS_noNT , fESC_b)
CCS_nofESC_noNT , CCS_nofESC_NT= diff_boundary(CCS_nofESC , common_NT_b)
NT_nofESC_noCCS , NT_nofESC_CCS = diff_boundary(NT_nofESC , CCS_b)

noCCS_NT_nofESC , noCCS_NT_fESC = diff_boundary(noCCS_NT , fESC_b)
noCCS_fESC_noNT , noCCS_fESC_NT = diff_boundary(noCCS_fESC , common_NT_b)
noNT_fESC_noCCS , noNT_fESC_CCS = diff_boundary(noNT_fESC , CCS_b)


Abc = len(CCS_b) - len(CCS_NT) - len(CCS_fESC) + len(CCS_NT_fESC)
aBc = len(common_NT_b) - len(CCS_NT) - len(NT_fESC) + len(CCS_NT_fESC)
ABc = len(CCS_NT) - len(CCS_NT_fESC)
abC = len(fESC_b) - len(CCS_fESC) - len(NT_fESC) + len(CCS_NT_fESC)
AbC = len(CCS_fESC) - len(CCS_NT_fESC)
aBC = len(NT_fESC) - len(CCS_NT_fESC)
ABC = len(CCS_NT_fESC)

print Abc , aBc , ABc , abC , AbC , aBC , ABC


Abc = len(CCS_noNT_nofESC)
aBc = len(NT_nofESC_noCCS)
ABc = len(CCS_nofESC_NT)
abC = len(noCCS_fESC_noNT)
AbC = len(CCS_noNT_fESC)
aBC = len(noCCS_fESC_NT)
ABC = len(CCS_NT_fESC)

print Abc , aBc , ABc , abC , AbC , aBC , ABC
 

boundaries = {'CCS_noNT_nofESC' : CCS_noNT_nofESC , 'CCS_noNT_fESC' : CCS_noNT_fESC , 'CCS_NT_nofESC' : CCS_nofESC_NT,
              'noCCS_NT_fESC' : noCCS_fESC_NT , 'noCCS_NT_nofESC' : NT_nofESC_noCCS , 'noCCS_noNT_fESC' : noCCS_fESC_noNT ,
              'CCS_NT_fESC' : CCS_NT_fESC}


##Insulation_score
ins_type = np.dtype({'names' : ['chr' , 'site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'] , 
                     'formats' : ['S8' , np.int , np.float , np.float , np.float , np.float]})
                     
CCS_noNT_nofESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\CCS_noNT_nofESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)
noCCS_NT_fESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\noCCS_NT_fESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)
noCCS_NT_nofESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\noCCS_NT_nofESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)
noCCS_noNT_fESC_ins = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Insulation_score\\noCCS_noNT_fESC_Insulation_score_40K.txt' , skiprows = 1 , dtype = ins_type)


CCS_noNT_nofESC = []
for i in CCS_noNT_nofESC_ins:
    if (i['CCS'] != 0) and (i['CCS'] / i['NT5'] > 1.5) and (i['CCS'] / i['NT6'] > 1.5) and (i['CCS'] / i['fESC'] > 1.5):
        CCS_noNT_nofESC.append((i['chr'] , i['site']))
CCS_noNT_nofESC = np.array(CCS_noNT_nofESC , dtype = boun_type)

noCCS_fESC_NT = []
for i in noCCS_NT_fESC_ins:
    if (i['NT5'] != 0) and (i['NT6'] != 0) and (i['fESC'] != 0) and (i['NT5'] / i['CCS'] > 1) and (i['NT6'] / i['CCS'] > 1) and (i['fESC'] / i['CCS'] > 1):
        noCCS_fESC_NT.append((i['chr'] , i['site']))
noCCS_fESC_NT = np.array(noCCS_fESC_NT , dtype = boun_type)


NT_nofESC_noCCS = []
for i in noCCS_NT_nofESC_ins:
    if (i['NT5'] != 0) and (i['NT6'] != 0) and (i['fESC'] != 0) and (i['NT5'] / i['CCS'] > 1) and (i['NT6'] / i['CCS'] > 1) and (i['fESC'] / i['CCS'] > 1):
        NT_nofESC_noCCS.append((i['chr'] , i['site']))
NT_nofESC_noCCS = np.array(NT_nofESC_noCCS , dtype = boun_type)
    
noCCS_fESC_noNT = []
for i in noCCS_noNT_fESC_ins:
    if (i['NT5'] != 0) and (i['NT6'] != 0) and (i['fESC'] != 0) and (i['NT5'] / i['CCS'] > 1) and (i['NT6'] / i['CCS'] > 1) and (i['fESC'] / i['CCS'] > 1):
        noCCS_fESC_noNT.append((i['chr'] , i['site']))
noCCS_fESC_noNT = np.array(noCCS_fESC_noNT , dtype = boun_type)



Abc = len(CCS_noNT_nofESC)
aBc = len(NT_nofESC_noCCS)
ABc = len(CCS_nofESC_NT)
abC = len(noCCS_fESC_noNT)
AbC = len(CCS_noNT_fESC)
aBC = len(noCCS_fESC_NT)
ABC = len(CCS_NT_fESC)

print Abc , aBc , ABc , abC , AbC , aBC , ABC

boundaries = {'CCS_noNT_nofESC' : CCS_noNT_nofESC , 'CCS_noNT_fESC' : CCS_noNT_fESC , 'CCS_NT_nofESC' : CCS_nofESC_NT,
              'noCCS_NT_fESC' : noCCS_fESC_NT , 'noCCS_NT_nofESC' : NT_nofESC_noCCS , 'noCCS_noNT_fESC' : noCCS_fESC_noNT ,
              'CCS_NT_fESC' : CCS_NT_fESC}



fig = plt.figure(figsize = (10, 10))
venn3(subsets=(Abc , aBc , ABc , abC , AbC , aBC , ABC), set_labels=('CCS', 'NT' , 'fESC'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S8_figs\\CCS_NT(common)fESC_Boundary_ChrX_Venn3.pdf')

for cl in boundaries:
    out = open('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_Venn3_classify\\Boundary_' + cl + '_chrX.txt' , 'w')
    for i in boundaries[cl]:
        out.writelines(i['chr'] + '\t' + str(i['site']) + '\n')
    out.close()



common2 , n1 , n2 = Common2_boundary(common_NT_b , fESC_b)

Ab = len(common_NT_b) - len(common2)
aB = len(fESC_b) - len(common2)
AB = len(common2)

print Ab , aB , AB

fig = plt.figure(figsize = (10, 10))
venn2(subsets = (Ab , aB , AB), set_labels = ('CCS', 'fESC'))

run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S8_figs\\common_NT_fESC_ChrX_Venn2.pdf')

    

##-----------------------------------Boundary_genes---------------------------------------------------
union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)



genes = {}

for cl in boundaries:
    genes[cl] = []
    for i in boundaries[cl]:
        chro = i['chr']
        start = i['site'] - 40000
        end = i['site'] + 40000 
        tmp_gene = union[union['chr'] == 'chr' + chro]        
        mask = (tmp_gene['gene_site'] >= start) & (tmp_gene['gene_site'] <= end)
        overlap = tmp_gene[mask]
        if overlap.size != 0 :
            for j in overlap:
                genes[cl].append(j)
    genes[cl] = np.array(genes[cl] , dtype = union_type)
                
for k,v in genes.items():
    print k , len(v)
    Write2genes('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\' + k + '_genes_ChrX.txt' , v)                


CCS_speci = genes['CCS_noNT_nofESC']
common3 = genes['CCS_NT_fESC']
ESC_speci = np.hstack([genes['noCCS_NT_fESC'] , genes['noCCS_NT_nofESC'] , genes['noCCS_noNT_fESC']])

fig1 = Box_plot([np.log2(CCS_speci['CCS'] + 1) , np.log2(CCS_speci['NT5'] + 1) , np.log2(CCS_speci['NT6'] + 1) , np.log2(CCS_speci['fESC']+1)])

fig2 = Box_plot([np.log2(common3['CCS'] + 1) , np.log2(common3['NT5'] + 1) , np.log2(common3['NT6'] + 1) , np.log2(common3['fESC'] + 1)] )

fig3 = Box_plot([np.log2(ESC_speci['CCS'] + 1), np.log2(ESC_speci['NT5'] + 1) , np.log2(ESC_speci['NT6'] +1), np.log2(ESC_speci['fESC'] + 1)])

run_Plot(fig1 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig5_figs\\CCS_specific_boundary_genes.pdf')
run_Plot(fig2 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig5_figs\\Common3_boundary_genes.pdf')
run_Plot(fig3 , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\Fig5_figs\\fESC_specific_boundary_genes.pdf')



##-----------------------------------Boundary_genes_Venn2----------------------------------------------

common2 = np.hstack((boundaries['noCCS_NT_fESC'] , boundaries['CCS_NT_fESC']))


genes = []
for i in common2:
    chro = i['chr']
    start = i['site'] - 40000
    end = i['site'] + 40000
    tmp = union[union['chr'] == 'chr' + chro]
    mask = (tmp['gene_site'] >= start) & (tmp['gene_site'] <= end)
    overlap = tmp[mask]
    if overlap.size != 0:
        for j in overlap:
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
        
percent_b = len(diff) / len(genes)


##-----------------------------------Boundary_RNA_RPKM---------------------------------------------------
CCS_RNA = Sig_To_1K('H:\\Workspace_New\\data\\RNA\\signal\\normalization\\bedgraph_1K\\CCS_RNA_1K.bedgraph')
NT5_RNA = Sig_To_1K('H:\\Workspace_New\\data\\RNA\\signal\\normalization\\bedgraph_1K\\NT5_RNA_1K.bedgraph')
NT6_RNA = Sig_To_1K('H:\\Workspace_New\\data\\RNA\\signal\\normalization\\bedgraph_1K\\NT6_RNA_1K.bedgraph')
fESC_RNA = Sig_To_1K('H:\\Workspace_New\\data\\RNA\\signal\\normalization\\bedgraph_1K\\fESC_RNA_1K.bedgraph')
  
for k in CCS_RNA.keys():
    print k , len(CCS_RNA[k]) , len(NT5_RNA[k]) , len(NT6_RNA[k]) , len(fESC_RNA[k]) 
    
union_RNA = []
union_RNA_type = np.dtype({'names':['chr' , 'start' , 'end' , 'CCS' , 'NT5' , 'NT6' , 'fESC'],
                           'formats':['S8' , np.int , np.int , np.float , np.float , np.float , np.float]})

al_RNA = {'CCS' : CCS_RNA , 'NT5' : NT5_RNA , 'NT6' : NT6_RNA , 'fESC' : fESC_RNA}
chrom = ['X']



for g in chrom:
    for i in range(len(CCS_RNA[g])):
        chro = g
        start = i * 1000
        end = (i + 1) * 1000
        ccs = al_RNA['CCS'][g][i]
        nt5 = al_RNA['NT5'][g][i]
        nt6 = al_RNA['NT6'][g][i]
        fesc = al_RNA['fESC'][g][i]
        union_RNA.append((chro , start , end , ccs , nt5 , nt6 , fesc))
    
    
union_RNA = np.array(union_RNA , dtype = union_RNA_type)    
    
    

RNA = {}

for cl in ['CCS_noNT_nofESC' , 'noCCS_NT_fESC' , 'noCCS_noNT_fESC' , 'noCCS_NT_nofESC' , 'CCS_NT_fESC']:
    RNA[cl] = []
    for i in boundaries[cl]:
        chro = i['chr']
        start = i['site'] - 40000
        end = i['site'] + 40000 
        ccs = []
        nt5 = []
        nt6 = []
        fesc = []
        tmp_RNA = union_RNA[union_RNA['chr'] == chro]        
        mask = (tmp_RNA['start'] >= start) & (tmp_RNA['start'] <= end)
        overlap = tmp_RNA[mask]
        if overlap.size != 0 :
            for j in overlap:
                ccs.append(j['CCS'])
                nt5.append(j['NT5'])
                nt6.append(j['NT6'])
                fesc.append(j['fESC'])

        if sum(ccs[:40]) >= sum(ccs[41:]):
            ccs = sum(ccs[:40]) / sum(ccs[41:])
        else:
            ccs = sum(ccs[41:]) / sum(ccs[:40])
            
        if sum(nt5[:40]) >= sum(nt5[41:]):
            nt5 = sum(nt5[:40]) / sum(nt5[41:])
        else:
            nt5 = sum(nt5[41:]) / sum(nt5[:40])
        
        if sum(nt6[:40]) >= sum(nt6[41:]):
            nt6 = sum(nt6[:40]) / sum(nt6[41:])
        else:
            nt6 = sum(nt6[41:]) / sum(nt6[:40])
            
        if sum(fesc[:40]) >= sum(fesc[41:]):
            fesc = sum(fesc[:40]) / sum(fesc[41:])
        else:
            fesc = sum(fesc[41:]) / sum(fesc[:40])
        RNA[cl].append((chro , start , end , ccs , nt5 , nt6 , fesc))
    RNA[cl] = np.array(RNA[cl] , dtype = union_RNA_type)
    
RNA_new = {}
for k , v in RNA.items():
    RNA_new[k] = []
    for i in v:
        tmp = i['CCS'] + i['NT5'] + i['NT6'] + i['fESC']
        if np.isnan(tmp) or np.isinf(tmp):
            pass
        else:
            RNA_new[k].append(i)
    RNA_new[k] = np.array(RNA_new[k] , dtype = union_RNA_type)
    
CCS_speci = RNA_new['CCS_noNT_nofESC']
common3 = RNA_new['CCS_NT_fESC']
ESC_speci = np.hstack([RNA_new['noCCS_NT_fESC'] , RNA_new['noCCS_NT_nofESC'] , RNA_new['noCCS_noNT_fESC']])

Box_plot([CCS_speci['CCS'] , CCS_speci['NT5'] , CCS_speci['NT6'] , CCS_speci['fESC']] )

Box_plot([common3['CCS'] , common3['NT5'] , common3['NT6'] , common3['fESC']] )

Box_plot([ESC_speci['CCS'] , ESC_speci['NT5'] , ESC_speci['NT6'] , ESC_speci['fESC']])             



##-----------------------------Venn3_NT5_NT6_fESC-------------------------------------

union_NT5_NT6_fESC = union_boundary([NT5_b , NT6_b , fESC_b])

NT5_NT6 , NT5_noNT6 , noNT5_NT6 = Common2_boundary(NT5_b , NT6_b)
NT5_fESC , NT5_nofESC , noNT5_fESC  = Common2_boundary(NT5_b , fESC_b)
NT6_fESC , NT6_nofESC , noNT6_fESC = Common2_boundary(NT6_b , fESC_b)
NT5_NT6_fESC = Common_boundary([NT5_b , NT6_b , fESC_b] , union_NT5_NT6_fESC)

NT5_noNT6_nofESC , NT5_noNT6_fESC = diff_boundary(NT5_noNT6 , fESC_b)
NT5_nofESC_noNT6 , NT5_nofESC_NT6 = diff_boundary(NT5_nofESC , NT6_b)
NT6_nofESC_noNT5 , NT6_nofESC_NT5 = diff_boundary(NT6_nofESC , NT5_b)

noNT5_NT6_nofESC , noNT5_NT6_fESC = diff_boundary(noNT5_NT6 , fESC_b)
noNT5_fESC_noNT6 , noNT5_fESC_NT6 = diff_boundary(noNT5_fESC , NT6_b)
noNT6_fESC_noNT5 , noNT6_fESC_NT5 = diff_boundary(noNT6_fESC , NT5_b)


Abc = len(NT5_b) - len(NT5_NT6) - len(NT5_fESC) + len(NT5_NT6_fESC)
aBc = len(NT6_b) - len(NT5_NT6) - len(NT6_fESC) + len(NT5_NT6_fESC)
ABc = len(NT5_NT6) - len(NT5_NT6_fESC)
abC = len(fESC_b) - len(NT5_fESC) - len(NT6_fESC) + len(NT5_NT6_fESC)
AbC = len(NT5_fESC) - len(NT5_NT6_fESC)
aBC = len(NT6_fESC) - len(NT5_NT6_fESC)
ABC = len(NT5_NT6_fESC)

print Abc , aBc , ABc , abC , AbC , aBC , ABC


Abc = len(NT5_noNT6_nofESC)
aBc = len(NT6_nofESC_noNT5)
ABc = len(NT6_nofESC_NT5)
abC = len(noNT5_fESC_noNT6)
AbC = len(NT5_noNT6_fESC)
aBC = len(noNT5_NT6_fESC)
ABC = len(NT5_NT6_fESC)

print Abc , aBc , ABc , abC , AbC , aBC , ABC

fig = plt.figure(figsize = (10, 10))
venn3(subsets=(Abc , aBc , ABc , abC , AbC , aBC , ABC), set_labels=('NT5', 'NT6' , 'fESC'))
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Plot\\S8_figs\\NT5_NT6_fESC_Boundary_ChrX_Venn3.pdf')


boundaries = {'NT5_noNT6_nofESC' : NT5_noNT6_nofESC , 'NT5_noNT6_fESC' : NT5_noNT6_fESC , 'NT5_NT6_nofESC' : NT6_nofESC_NT5,
              'noNT5_NT6_fESC' : noNT5_NT6_fESC , 'noNT5_NT6_nofESC' : NT6_nofESC_noNT5 , 'noNT5_noNT6_fESC' : noNT5_fESC_noNT6 ,
              'NT5_NT6_fESC' : NT5_NT6_fESC}

union_gene = get_union_gene('H:\\Workspace_New\\data\\RNA\\gene_expression\\all_gene_expression.txt')

union = []
for i in union_gene:
    if (i['CCS'] > 1) or (i['NT5'] > 1) or (i['NT6'] > 1) or (i['fESC'] > 1):
        gene_site = int((i['start'] + i['end']) / 2)
        union.append((i['gene_name'] , i['chr'] , gene_site , i['CCS'] , i['NT5'] , i['NT6'] , i['fESC']))
union = np.array(union , dtype = union_type)



genes = {}

for cl in boundaries:
    genes[cl] = []
    for i in boundaries[cl]:
        chro = i['chr']
        start = i['site'] - 40000
        end = i['site'] + 40000 
        tmp_gene = union[union['chr'] == 'chr' + chro]        
        mask = (tmp_gene['gene_site'] >= start) & (tmp_gene['gene_site'] <= end)
        overlap = tmp_gene[mask]
        if overlap.size != 0 :
            for j in overlap:
                genes[cl].append(j)
    genes[cl] = np.array(genes[cl] , dtype = union_type)
                
for k,v in genes.items():
    print k , len(v)
    Write2genes('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_ChrX_Venn3(NT5_NT6_fESC)\\' + k + '_genes_ChrX.txt' , v)                



