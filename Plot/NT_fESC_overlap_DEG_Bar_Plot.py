# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:57:01 2020

@author: han-luo
"""
from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt



gene_type = np.dtype({'names': ['gene_name' , 'chr' , 'gene_site' , 'CCS' , 'NT5' , 'NT6' , 'fESC'] , 
                     'formats': ['S64' , 'S8' , np.int , np.float , np.float , np.float , np.float]})
                     
##Boundary                     
data1 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\CCS_NT_nofESC_genes.txt' , dtype = gene_type , skiprows = 1)
data2 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\noCCS_NT_nofESC_genes.txt' , dtype = gene_type , skiprows = 1)
data3 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\CCS_noNT_fESC_genes.txt' , dtype = gene_type , skiprows = 1)
data4 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\noCCS_noNT_fESC_genes.txt' , dtype = gene_type , skiprows = 1)

diff1 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\CCS_NT_nofESC_diff_genes.txt' , dtype = gene_type , skiprows = 1)
diff2 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\noCCS_NT_nofESC_diff_genes.txt' , dtype = gene_type , skiprows = 1)
diff3 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\CCS_noNT_fESC_diff_genes.txt' , dtype = gene_type , skiprows = 1)
diff4 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\HiTAD\\respective\\stable_600K\\bottom_domain\\Boundary_Venn3\\boundary_genes\\noCCS_noNT_fESC_diff_genes.txt' , dtype = gene_type , skiprows = 1)


##Loop
data5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\CCS_NT_nofESC_genes.txt' , dtype = gene_type , skiprows = 1)
data6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\noCCS_NT_nofESC_genes.txt' , dtype = gene_type , skiprows = 1)
data7 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\CCS_noNT_fESC_genes.txt' , dtype = gene_type , skiprows = 1)
data8 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\noCCS_noNT_fESC_genes.txt' , dtype = gene_type , skiprows = 1)

diff5 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\CCS_NT_nofESC_diff_genes.txt' , dtype = gene_type , skiprows = 1)
diff6 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\noCCS_NT_nofESC_diff_genes.txt' , dtype = gene_type , skiprows = 1)
diff7 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\CCS_noNT_fESC_diff_genes.txt' , dtype = gene_type , skiprows = 1)
diff8 = np.loadtxt('H:\\Workspace_New\\data\\HiC\\Loop\\Raw_20K_0.05\\union_Loops\\Loop_Venn3_genes\\noCCS_noNT_fESC_diff_genes.txt' , dtype = gene_type , skiprows = 1)


percent_b_NT = (len(diff1) + len(diff2)) / (len(data1) + len(data2))
percent_b_fESC = (len(diff3) + len(diff4)) / (len(data3) + len(data4))

percent_l_NT = (len(diff5) + len(diff6)) / (len(data5) + len(data6))
percent_l_fESC = (len(diff7) + len(diff8)) / (len(data7) + len(data8))




def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2.- 0.2, 1.03*height, '%s' % np.round(height , 4))


name_list = ['Boundary_NT', 'Boundary_fESC' , 'Loop_NT' , 'Loop_fESC']
num_list = [percent_b_NT, percent_b_fESC , percent_l_NT, percent_l_fESC]


autolabel(plt.bar(range(len(num_list)), num_list, color='plum', tick_label=name_list))
plt.show()