# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 15:29:06 2019

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
    
def Bar_plot(x , y , p_list , color):
    left, bottom, width, height = 0.45, 0.1, 0.5, 0.8
    size_axes = [left, bottom, width, height]
    fig = plt.figure(figsize = (12, 12))
    ax = fig.add_axes(size_axes)
    numbers = [len(i) for i in p_list]
    colors = []
    for i in range(len(numbers)):
        colors.extend([color[i] for j in range(numbers[i] + 1)])
    ax.barh(x , y , color = colors , label = "cluster3")
    yticks = [i for i in x]
    labels = []
    for i in range(len(numbers)):
        for j in p_list[i]:
            labels.append(j)
        labels.append('')
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
goFil_ESC = 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\3 go.xlsx'

pathway_cluster2 = get_datas(goFil_CCS)
pathway_cluster3 = get_datas(goFil_ESC)
#CCS
colors_CCS = ['yellow' , 'gray' , 'skyblue' , 'lightgreen' ,'coral']
pathway_CCS1 = ['transport' , 'cell adhesion' , 'cell differentiation' ]
pathway_CCS2 = ['lipid metabolic process' , 'fatty acid metabolic process' ,'regulation of cholesterol metabolic process' , 
                'sphingolipid metabolic process' , 'retinal metabolic process' , 'metabolic process' , 'arachidonic acid metabolic process' , 
                'xenobiotic metabolic process' , 'prostaglandin metabolic process' , 'retinoic acid metabolic process']
pathway_CCS3 = ['positive regulation of osteoblast differentiation' , 'smooth muscle cell differentiation' , 
                'chondrocyte differentiation' , 'muscle cell differentiation' , 'negative regulation of fat cell differentiation' ,
                'positive regulation of melanocyte differentiation' , 'epithelial cell differentiation' ,'glial cell differentiation' ,
                'positive regulation of chondrocyte differentiation' , 'positive regulation of natural killer cell differentiation' ,
                'regulation of fat cell differentiation' , 'regulation of epithelial cell differentiation' , 'positive regulation of fat cell differentiation' ,
                'inner ear receptor cell differentiation' ]
pathway_CCS4 = ['neuron projection development' , 'chondrocyte development' , 'mammary gland alveolus development' ,
                'multicellular organism development' , 'male gonad development' , 'cartilage development' , 'nervous system development' ,
                'uterus development' , 'lung alveolus development' , 'heart development' , 'positive regulation of cartilage development' ,
                'inner ear development' , 'diaphragm development' , 'skeletal muscle tissue development' , 'exocrine pancreas development' ,
                'cell growth involved in cardiac muscle cell development' , 'ovarian follicle development' , 'positive regulation of neuron projection development' ,
                'cardiovascular system development' , 'skeletal system development' , 'ventricular septum development' , 'negative regulation of female gonad development' ,
                'development of secondary male sexual characteristics' , 'aorta development' , 'ureteric bud development']
pathway_CCS5 = ['positive regulation of MAPK cascade' , 'activation of MAPK activity' , 'activation of MAPKK activity']

CCS_list = [pathway_CCS1 , pathway_CCS2 , pathway_CCS3 , pathway_CCS4 , pathway_CCS5]
x_CCS , y_CCS = get_selected_p_value(CCS_list , pathway_cluster2)
fig = Bar_plot(x_CCS , y_CCS , CCS_list , colors_CCS)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\cluster2_go.pdf')

#ESC
colors_ESC = ['yellow' , 'gray' , 'skyblue' , 'lightgreen' ,'coral' , 'navajowhite']
pathway_ESC1 = ['multicellular organism development' , 'DNA replication']
pathway_ESC2 = ['embryonic digestive tract morphogenesis' , 'embryonic hindlimb morphogenesis' , 'in utero embryonic development' ,
                'embryonic skeletal system morphogenesis' , 'embryonic limb morphogenesis' , 'embryonic cranial skeleton morphogenesis' ,
                'embryonic forelimb morphogenesis' , 'embryonic skeletal joint morphogenesis' , 'embryonic digit morphogenesis' ,
                'embryonic pattern specification' , 'embryonic camera-type eye morphogenesis' , 'embryonic placenta development']
pathway_ESC3 = ['cell cycle' , 'cell division' , 'mitotic nuclear division' , 'meiotic cell cycle' , 'chromosome segregation' ,
                'mitotic cell cycle' , 'reciprocal meiotic recombination' , 'positive regulation of mitotic cell cycle' ,
                'regulation of mitotic cell cycle' , 'meiotic nuclear division' , 'positive regulation of cell cycle' ,
                'regulation of G1/S transition of mitotic cell cycle' , 'positive regulation of G2/M transition of mitotic cell cycle' ,
                'G1/S transition of mitotic cell cycle' , 'meiotic sister chromatid cohesion, centromeric']
pathway_ESC4 = ['cellular response to DNA damage stimulus' , 'DNA repair' , 'intrinsic apoptotic signaling pathway in response to DNA damage' ,
                'intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator' , 'negative regulation of DNA damage response, signal transduction by p53 class mediator']
pathway_ESC5 = ['somatic stem cell population maintenance' , 'neuronal stem cell population maintenance' , 'stem cell population maintenance' ,
                'positive regulation of stem cell proliferation' , 'germ-line stem cell population maintenance']
pathway_ESC6 = ['telomere maintenance' , 'telomere maintenance via recombination']

ESC_list = [pathway_ESC1 , pathway_ESC2 , pathway_ESC3 , pathway_ESC4 , pathway_ESC5 , pathway_ESC6]
x_ESC , y_ESC = get_selected_p_value(ESC_list , pathway_cluster3)
fig = Bar_plot(x_ESC , y_ESC , ESC_list , colors_ESC)
run_Plot(fig , 'D:\\ntESC_3Dreprogramming\\Figures\\Fig1\\new\\cluster3_go.pdf')

