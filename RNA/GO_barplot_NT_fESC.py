# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 17:47:33 2019

@author: han-luo
"""

import xlrd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_datas(goFil):
    d_type = ({'names':['p_name'  , 'p_value'],
               'formats':['S128' , np.float]})
    goData = xlrd.open_workbook(goFil)
    sheet = goData.sheets()[0]
    pathway = sheet.col_values(1)
    p_value = sheet.col_values(4)
    goDatas = []
    for i in range(len(pathway)):
        if '~' not in pathway[i]:
            continue
        p_name = pathway[i].strip().split('~')[1]
        p_v = -np.log10(p_value[i])
        goDatas.append((p_name , p_v) )
    goDatas = np.array(goDatas , dtype = d_type)
    return goDatas
    
def get_selected_p_value(p_list , p_classify):
    y = []
    for i in p_list:
        for j in i:
            p_v = p_classify[p_classify['p_name'] == j]['p_value'][0]
            y.append(p_v)
        y.append(0)
    x = range(len(y))
    return x , y
    
def Bar_plot(x , y):
    left, bottom, width, height = 0.3, 0.1, 0.6, 0.8
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    ax.barh(x , y , color = 'deepskyblue' , label = "Diff expression")
    yticks = [i for i in x]
    labels = pathway_CCS1
    ax.set_yticks(yticks)
    ax.set_yticklabels(labels,fontsize = 10)
    ax.set_xlabel('-log10(p value)' , fontsize = 20 )
    ax.legend(loc = 'upper right')
    return fig
    
def run_Plot(fig , OutFile):
    pp = PdfPages(OutFile)
    pp.savefig(fig)
    pp.close()
    

    
#GO_Files
goFil_CCS = 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\2 go.xlsx'

pathway_cluster2 = get_datas(goFil_CCS)
#CCS
colors_CCS = ['skyblue']
pathway_CCS1 = ['response to mechanical stimulus' , 'response to drug' , 'female pregnancy' , 'ossification' , 
                'cartilage development' , 'skeletal system development' , 'angiogenesis' , 'osteoblast differentiation' , 
                'response to progesterone' , 'response to estrogen' , 'response to estradiol']
y_CCS = np.array([6.51600914688057e-13 , 8.57302408188616E-07 , 0.0026174539249271 , 2.59373236794074E-08 , 0.000025224645903225 , 
                 0.000112726229411862 , 1.06889581890206E-06 , 0.00125325713372839 , 5.77700231711485E-08, 0.000190810955378653 , 0.000756496349797145])
y_CCS = -np.log10(y_CCS)
x_CCS = np.arange(len(y_CCS))

fig = Bar_plot(x_CCS , y_CCS)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\Diff_gene_expression_go.pdf')

