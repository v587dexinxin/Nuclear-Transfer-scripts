# -*- coding: utf-8 -*-
"""
Created on Wed May 23 20:52:50 2018

@author: han-luo
"""

from __future__ import division
from itertools import islice
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os

"""
    This is an Only HiC-ploting scripts, We will split the HiC Contact Matrix by 
    chromosome.For each chromosome will will split the whole Contact Map into chunk
    maps,That means only part of contact matrix consists of not far interactions.
    This map feature can be adjusted by parameters.
    
    Notice:
        Parameters
        ----------
        
"""


def Cmap_setting(**kwargs):
    
    """
        Color bar setting. hexadecimal notation(#16) must be input.
        
        Parameters
        ----------
        start_Color : str
            start_Color (default : #FFFFFF)        
        
        end_Color : str
            end_Color (default : #CD0000)
            
    """
    start_Color = kwargs.get('start_Color','#FFFFFF')
    end_Color = kwargs.get('end_Color','#CD0000')
    
    return LinearSegmentedColormap.from_list('interactions',[start_Color,end_Color])


def getmatrix(inter,l_bin,r_bin,data_type):
    
    """
        Get 2D Interaction Matrix.
        Parameters
        ----------
        inter : numpy
            whole data for a chromosome
        
        l_bin : int
            Left bin
        
        r_bin : int
            Right bin
        
        data_type : str 
            A : 3 tuples Matrix for each chromosome in npz dict.
            B : whole 2D Matrix for each chromosome in npz dict.
            
    """
    if data_type == 'A':
        
        inter_matrix = np.zeros((r_bin - l_bin, r_bin - l_bin),dtype = float )
        #Extract the regional data
        mask = (inter['bin1'] >= l_bin) & (inter['bin1'] < r_bin) & \
               (inter['bin2'] >= l_bin) & (inter['bin2'] < r_bin)
        inter_extract = inter[mask]
    
        #Fill the matrix:
        for i in inter_extract:
            #Off-diagnoal parts
            if i['bin1'] != i['bin2']:
                inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
                inter_matrix[i['bin2'] - l_bin][i['bin1'] - l_bin] += i['IF']
            else:
                #Diagonal part
                inter_matrix[i['bin1'] - l_bin][i['bin2'] - l_bin] += i['IF']
    else:
        inter_matrix = inter[l_bin:r_bin,l_bin:r_bin]
    return inter_matrix


  
def Matplotlib_setting():
    """
        run Matplotlib settings.
    """
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'

def caxis_H(ax):
    """
        Axis Control for HeatMaps
    """
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = 10, length = 5, pad = 7)

def properU(pos):
    """
    Express a genomic position in a proper unit (KB, MB, or both).
    
    """
    i_part = int(pos) // 1000000 # Integer Part
    d_part = (int(pos) % 1000000) // 1000 # Decimal Part
    
    if (i_part > 0) and (d_part > 0):
        return ''.join([str(i_part), 'M', str(d_part), 'K'])
    elif (i_part == 0):
        return ''.join([str(d_part), 'K'])
    else:
        return ''.join([str(i_part), 'M'])

def Load_HiC_Data(npz):
    """
        return HiC source.Input must be npz
    """
    return np.load(npz)

def Load_loops(loop_fil):
    """
    """
    loop_type = np.dtype({'names':['chr','start','end'],
                          'formats':['S4',np.int,np.int]})
    
    loops = np.loadtxt(loop_fil, dtype = loop_type, skiprows = 1, usecols = [0,1,2])
    
    return loops


def Load_5K_loops(loop_fil):
    """
    """
    loop_type = np.dtype({'names':['chr','start','end'],
                          'formats':['S4',np.int,np.int]})
    
    loops = np.loadtxt(loop_fil, dtype = loop_type, skiprows = 1, usecols = [0,1,4])
    
    return loops


def Load_Addtional_Loop(loop_fil):
    """
    """
    loop_type = np.dtype({'names':['chr','start','end','P_V','Sig'],
                          'formats':['S4',np.int,np.int,np.float,np.float]})
    loops = np.loadtxt(loop_fil,dtype = loop_type)
    
    return loops

def Load_loops_qvalues(loop_fil,chrom):
    """
    """
    loop_type = np.dtype({'names':['chr','start','end','qvalue'],
                          'formats':['S4',np.int,np.int,np.float]}) 

    loops = np.loadtxt(loop_fil, dtype = loop_type, skiprows = 1, usecols = [0,1,2,-1])
    
    mask = (loops['qvalue'] <= 0.05) & (loops['chr'] == chrom)
    loops = loops[mask]
    
    return loops

def Load_cluster_loops(loop_fil):
    """
    """
    loop_type = np.dtype({'names':['chr','start','end','wq'],
                          'formats':['S4',np.int,np.int,np.float]})
    
    loops = np.loadtxt(loop_fil, dtype = loop_type, skiprows = 1, usecols = [0,1,2,3])
    
    return loops

def Load_Allelic_Loops(loop_fil):
    """
    """
    loop_type = np.dtype({'names':['chr','start','end'],
                          'formats':['S4',np.int,np.int]})
    
    loops = np.loadtxt(loop_fil,dtype = loop_type,skiprows = 1,usecols = [0,1,2])
    
    return loops





def UpdateDI(DI):
    """
    """
    New_DI = []
    New_index = []

    for index in range(len(DI) - 1):
        if DI[index] * DI[index + 1] < 0:
            New_DI.append(DI[index])
            New_DI.append(0)
            New_index.append(index)
            New_index.append(index + 0.5)
        else:
            New_DI.append(DI[index])
            New_index.append(index)
    
    return np.array(New_index), np.array(New_DI)

def plotting_settings(**kwargs):
    """
        Plotting settings.
        Parameters
        ----------
        figsize : tuple
            about figure size (default : (12,14))
        
        width : float
            about figure width control (default : 0.8)
            
        HB : float
            about figure
        
        ResHiC : int
            HiC data Resolution
        
        interval : int
            Step Length.
            
    """
    figsize = kwargs.get('figsize',(12,14))
    width = kwargs.get('width',0.8)
    HB = kwargs.get('HB',0.1)
    ResHiC = kwargs.get('ResHiC',40000)
    interval = kwargs.get('interval',200)
    chroms = kwargs.get('chroms',['1','2','3','4','5','6','7','8','9','10',
                                  '11','12','13','14','15','16','17','18','19','X'])
    return {'figsize' : figsize,
            'width' : width,
            'HB' : HB,
            'ResHiC' : ResHiC,
            'interval' : interval,
            'chroms' : chroms}

def Plotting(prefix,Lib,res,interval,data_type = 'A'):
    """
        Plotting ...
        Parameters
        ----------
        prefix : str
            prefix of Output pdfs.
        
        npz : str
            HiC data.(npz)
        
        data_type : str
            A : 3 tuples Matrix for each chromosome in npz dict.
            B : whole 2D Matrix for each chromosome in npz dict.
    """
    
    #-------settings-----------
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    plot_sets = plotting_settings(ResHiC = res,interval = interval)
#    chrom = plot_sets['chroms']
    ResHiC = plot_sets['ResHiC']
    interval = plot_sets['interval']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]

    #------DataLoading---------
    #HiC_data = Load_HiC_Data(npz)
    
    #------Drawing---------
    chrom = Lib.keys()
    #chrom = ['1']
    for i in chrom:
        inter  = Lib[i]
        pp = PdfPages(prefix + '-' + str(i) + '.pdf')
        startHiC = 0
        if data_type == 'A':
            Len = inter['bin1'].max() // interval
        else:
            Len = inter.shape[0] // interval
        for idx in range(Len):
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left,HB,width,HH])
            EndHiC = startHiC + interval
            
            Matrix = getmatrix(inter,startHiC,EndHiC,data_type)
            nonzero = Matrix[np.nonzero(Matrix)]
            if nonzero.size <= 100:
                plt.close(fig)
                startHiC = EndHiC
                continue
            
            vmax = np.percentile(nonzero,95)
            
            sc = ax.imshow(Matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, interval, 0, interval), vmax = vmax, origin = 'lower')
            
            ticks = list(np.linspace(0,interval,5).astype(int))
            pos = [((startHiC + t) * ResHiC) for t in ticks]
            labels = [properU(p) for p in pos]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            
            ax = fig.add_axes([Left + width + 0.02, HB, 0.01, HH])
            fig.colorbar(sc,cax = ax)
            
            startHiC = EndHiC
            pp.savefig(fig)
            plt.close(fig)
        pp.close()


def Plotting_DI_TADs(prefix,Lib,Sigs,TADs,res,interval,data_type = 'A'):
    """
        Plotting ...
        Parameters
        ----------
        prefix : str
            prefix of Output pdfs.
        
        npz : str
            HiC data.(npz)
        
        data_type : str
            A : 3 tuples Matrix for each chromosome in npz dict.
            B : whole 2D Matrix for each chromosome in npz dict.
    """
    
    #-------settings-----------
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    plot_sets = plotting_settings(ResHiC = res,interval = interval)
#    chrom = plot_sets['chroms']
    ResHiC = plot_sets['ResHiC']
    interval = plot_sets['interval']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]

    #------DataLoading---------
    #HiC_data = Load_HiC_Data(npz)
    
    #------Drawing---------
    chrom = Lib.keys()
    #chrom = ['1']
    for i in chrom:
        inter  = Lib[i]
        Sigs_chr = Sigs[i]
        TADs_chr = TADs[i]
        pp = PdfPages(prefix + '-' + str(i) + '.pdf')
        startHiC = 0
        if data_type == 'A':
            Len = inter['bin1'].max() // interval
        else:
            Len = inter.shape[0] // interval
        for idx in range(Len):
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left,HB,width,HH])
            EndHiC = startHiC + interval
            
            Matrix = getmatrix(inter,startHiC,EndHiC,data_type)
            tmp_Sigs = Sigs_chr[startHiC:EndHiC]
            mask = ((TADs_chr['start'] > startHiC * ResHiC) & (TADs_chr['start'] < EndHiC * ResHiC)) |\
                   ((TADs_chr['end'] > startHiC * ResHiC) & (TADs_chr['end'] < EndHiC * ResHiC))
            tads = TADs_chr[mask]
            nonzero = Matrix[np.nonzero(Matrix)]
            if nonzero.size <= 100:
                plt.close(fig)
                startHiC = EndHiC
                continue
            
            vmax = np.percentile(nonzero,95)
            
            sc = ax.imshow(Matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, interval, 0, interval), vmax = vmax, origin = 'lower')
            
            for tad in tads:
                s_t = tad['start'] // res - startHiC
                s_e = tad['end'] // res - startHiC
                ax.plot([s_t,s_e],[s_t,s_t],color = 'black',ls = '--',lw = 1.0)
                ax.plot([s_e,s_e],[s_t,s_e],color = 'black',ls = '--',lw = 1.0)
            
            ax.set_xlim(0,(EndHiC-startHiC))
            ax.set_ylim(0,(EndHiC-startHiC))            
            
            ticks = list(np.linspace(0,interval,5).astype(int))
            pos = [((startHiC + t) * ResHiC) for t in ticks]
            labels = [properU(p) for p in pos]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            
            ax2 = fig.add_axes([Left,HB+HH,width,0.1])
            for spine in ['right','top','left']:
                ax2.spines[spine].set_visible(False)            
            index,sigs = UpdateDI(tmp_Sigs)
            ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
            ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
            ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                            right = False, labelbottom =False, labeltop = False,
                            labelleft = False, labelright = False)
            ax2.set_xlim(0,(EndHiC-startHiC))            
            startHiC = EndHiC
            pp.savefig(fig)
            plt.close(fig)
        pp.close()







def Plotting_SubArea(out,Lib,chro,start,end,res):
    """
        Plot Loop area.
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    Matrix = Lib[chro]
    s = start // res
    e = end // res
#    s1_s = s1 // res - s
#    e1_s = e1 // res - s
    M = Matrix[s:e,s:e]

    nonzero = M[np.nonzero(M)]
    vmax = np.percentile(nonzero,90)
#    print vmax
    fig,ax = plt.subplots(1,figsize = (10,9))
    ax.imshow(M.T,cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                  extent = (0,M.shape[0],0,M.shape[0]),vmax = vmax, origin = 'lower')
#    ax.plot([s1_s,e1_s],[s1_s,s1_s],color = 'black',lw = 0.5)
#    ax.plot([s1_s,e1_s],[e1_s,e1_s],color = 'black',lw = 0.5)
#    ax.plot([s1_s,s1_s],[s1_s,e1_s],color = 'black',lw = 0.5)
#    ax.plot([e1_s,e1_s],[s1_s,e1_s],color = 'black',lw = 0.5)
    
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    
    
    plt.savefig(out)
    plt.close()

def Plotting_SubArea_With_Sigs(out,Lib,Sigs,chro,start,end,res,s1,e1):
    """
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    Matrix = Lib[chro]
    s = start // res
    e = end // res
    s1_s = s1 // res - s
    e1_s = e1 // res - s
    M = Matrix[s:e,s:e]
    NewSigs = Sigs[s:e]
    nonzero = M[np.nonzero(M)]
    vmax = np.percentile(nonzero,96.5)
    size = (10,9)
    width = 0.618;Left = (1-width) / 2
    HB = 0.1;HH = width * size[0] / size[1]
    SB = HB + HH
    
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left,HB,width,HH])
    ax.imshow(M,cmap = my_cmap, aspect = 'auto',interpolation = 'none',
              extent = (0,M.shape[0],0,M.shape[0]),vmax = vmax, origin = 'lower')

    ax.plot([s1_s,e1_s],[s1_s,s1_s],color = 'black',lw = 0.5)
    ax.plot([s1_s,e1_s],[e1_s,e1_s],color = 'black',lw = 0.5)
    ax.plot([s1_s,s1_s],[s1_s,e1_s],color = 'black',lw = 0.5)
    ax.plot([e1_s,e1_s],[s1_s,e1_s],color = 'black',lw = 0.5)
    
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    ax2 = fig.add_axes([Left,SB,width,0.1])
    for spine in ['right','top','left']:
        ax2.spines[spine].set_visible(False)
    
    index,sigs = UpdateDI(NewSigs)
    ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
    ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
    ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    ax2.set_xlim(0,len(NewSigs))

    plt.savefig(out)
    plt.close(fig)

def Plotting_SubArea_With_Sigs_TADs(out,Lib,Sigs,TADs,chro,start,end,res):
    """
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    Matrix = Lib[chro]
    tads = TADs[chro]
    s = start // res
    e = end // res
    mask = ((tads['end'] > start) & (tads['end'] < end)) | \
            ((tads['start'] > start) & (tads['start'] < end))
    
#    mask = (tads['start'] > start) & (tads['end'] < end)
    sub_tads = tads[mask]
    print sub_tads
        
#    s1_s = s1 // res - s
#    e1_s = e1 // res - s
    M = Matrix[s:e,s:e]
    NewSigs = Sigs[s:e]
    nonzero = M[np.nonzero(M)]
    vmax = np.percentile(nonzero,95)
    size = (10,9)
    width = 0.618;Left = (1-width) / 2
    HB = 0.1;HH = width * size[0] / size[1]
    SB = HB + HH
    
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left,HB,width,HH])
    ax.imshow(M,cmap = my_cmap, aspect = 'auto',interpolation = 'none',
              extent = (0,M.shape[0],0,M.shape[0]),vmax = vmax, origin = 'lower')

    for tad in sub_tads:
        s_t = tad['start'] // res - s
        s_e = tad['end'] // res - s    
        ax.plot([s_t,s_e],[s_t,s_t],color = 'black', ls = '--',lw = 1.75)
        ax.plot([s_e,s_e],[s_t,s_e],color = 'black', ls = '--',lw = 1.75)
    
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    ax.set_xlim(0,(e-s))
    ax.set_ylim(0,(e-s))
    
    ax2 = fig.add_axes([Left,SB,width,0.1])
    for spine in ['right','top','left']:
        ax2.spines[spine].set_visible(False)
    
    index,sigs = UpdateDI(NewSigs)
    ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
    ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
    ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    ax2.set_xlim(0,(e-s))

    plt.savefig(out)
    plt.close(fig)




def Plotting_SubArea_With_Sigs_New(out,Lib,Sigs,chro,start,end,res,ls,le):
    """
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    Matrix = Lib[chro]
    s = start // res
    e = end // res
    M = Matrix[s:e,s:e]
    NewSigs = Sigs[s*40:e*40]
    ls_ = ls // 1000 - start // 1000
    le_ = le // 1000 - start // 1000
    NewSigs[:ls_-0] = 0
    NewSigs[ls_+80:le_-40] = 0
    NewSigs[le_+80:] = 0
    nonzero = M[np.nonzero(M)]
    vmax = np.percentile(nonzero,85)
    print vmax
    size = (10,9)
    width = 0.618;Left = (1-width) / 2
    HB = 0.1;HH = width * size[0] / size[1]
    SB = HB + HH
    
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left,HB,width,HH])
    ax.imshow(M,cmap = my_cmap, aspect = 'auto',interpolation = 'none',
              extent = (0,M.shape[0],0,M.shape[0]),vmax = 800, origin = 'lower')
    
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    ax2 = fig.add_axes([Left,SB,width,0.1])
    for spine in ['right','top','left']:
        ax2.spines[spine].set_visible(False)
    
    ax2.fill_between(np.arange(len(NewSigs)),NewSigs,color = 'b')
#    index,sigs = UpdateDI(NewSigs)
#    ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
#    ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
    ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    ax2.set_xlim(0,len(NewSigs))
    ax2.set_ylim(0,NewSigs.max()*0.7)
    plt.savefig(out)
    plt.close(fig)

def Plotting_SubArea_With_Sigs_New2(out,Lib,chro,start,end,res,loop):
    """
    """
    pp = PdfPages(out)
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    Matrix = Lib[chro]
    #title = '-'.join(map(str,loop[:4]))
    s = start // res
    e = end // res
    M = Matrix[s:e,s:e]
    ls = loop[1]
    le = loop[2]
    ls_ = ls // 40000
    le_ = le // 40000
    nonzero = M[np.nonzero(M)]
    vmax = np.percentile(nonzero,90)
    size = (10,9)
    width = 0.618;Left = (1-width) / 2
    HB = 0.1;HH = width * size[0] / size[1]
    SB = HB + HH
    
    fig = plt.figure(figsize = size)
    ax = fig.add_axes([Left,HB,width,HH])
    ax.imshow(M,cmap = my_cmap, aspect = 'auto',interpolation = 'none',
              extent = (0,M.shape[0],0,M.shape[0]),vmax = vmax, origin = 'lower')
    #ax.scatter(le // 40000 - s + 0.5,ls // 40000 - s + 0.5,color = '',edgecolors = 'b',s = 20)
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    ax2 = fig.add_axes([Left,SB,width,0.1])
    
    for spine in ['right','top','left']:
        ax2.spines[spine].set_visible(False)
    
    ax2.bar([ls_ - s  + 0.5,le_ - s + 0.5],[loop[-2],loop[-1]],color = 'b',width = 0.15)
#    index,sigs = UpdateDI(NewSigs)
#    ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
#    ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
    ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    ax2.set_xlim(0,(e-s))
    ax2.set_ylim(0,max(loop[-1],loop[-2])*0.7)
    #fig.suptitle(title)
    pp.savefig(fig)
    plt.close(fig)
    pp.close()

def Plotting_Whole(Loops,OutPath,Lib,Sigs,Allelic,dis = 2000000):
    """
    """
    pp = PdfPages(os.path.join(OutPath,'Allelic_Calling_'+Allelic+'.pdf'))
    for lp in Loops:
        chro = lp[0]
        start = lp[1]
        end = lp[2]
        
        tmp_dis = end - start
        if tmp_dis > dis:
            start_d = start - 400000
            end_d = end + 400000
        else:
            d  = (dis - tmp_dis) // 2
            start_d = start - d
            end_d = end + d
        
        Plotting_SubArea_With_Sigs_New2(pp,Lib,chro,start_d,end_d,40000,lp)
    
    pp.close()




def Plotting_SubArea_With_Sigs_Alleic_New(M_out,P_out,M_Lib,P_Lib,M_Sigs,P_Sigs,
                                   chro,start,end,res,ls,le):
    """
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    M_Matrix = M_Lib[chro]
    P_Matrix = P_Lib[chro]
    s = start // res
    e = end // res
    M_M = M_Matrix[s:e,s:e]
    P_P = P_Matrix[s:e,s:e]
    M_NewSigs = M_Sigs[s*40:e*40]
    P_NewSigs = P_Sigs[s*40:e*40]
    ls_ = ls // 1000 - start // 1000
    le_ = le // 1000 - start // 1000
    M_NewSigs[:ls_-40] = 0
    M_NewSigs[ls_+80:le_-40] = 0
    M_NewSigs[le_+40:] = 0
    P_NewSigs[:ls_-40] = 0
    P_NewSigs[ls_+80:le_-40] = 0
    P_NewSigs[le_+40:] = 0
    
    nonzero = M_M[np.nonzero(M_M)]
    #vmax = np.percentile(nonzero,90)
    size = (10,9)
    width = 0.618;Left = (1-width) / 2
    HB = 0.1;HH = width * size[0] / size[1]
    SB = HB + HH
    
    fig1 = plt.figure(figsize = size)
    ax = fig1.add_axes([Left,HB,width,HH])
    ax.imshow(M_M,cmap = my_cmap, aspect = 'auto',interpolation = 'none',
              extent = (0,M_M.shape[0],0,M_M.shape[0]),vmax = 44, origin = 'lower')
    
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    ax2 = fig1.add_axes([Left,SB,width,0.1])
    for spine in ['right','top','left']:
        ax2.spines[spine].set_visible(False)
    
    ax2.fill_between(np.arange(len(M_NewSigs)),M_NewSigs,color = 'b')
#    index,sigs = UpdateDI(NewSigs)
#    ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
#    ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
    ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    ax2.set_xlim(0,len(M_NewSigs))
    ax2.set_ylim(0,max(M_NewSigs.max(),P_NewSigs.max()))
    plt.savefig(M_out)
    plt.close(fig1)


    fig2 = plt.figure(figsize = size)
    ax = fig2.add_axes([Left,HB,width,HH])
    ax.imshow(P_P,cmap = my_cmap, aspect = 'auto',interpolation = 'none',
              extent = (0,P_P.shape[0],0,P_P.shape[0]),vmax = 44, origin = 'lower')
    
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    ax2 = fig2.add_axes([Left,SB,width,0.1])
    for spine in ['right','top','left']:
        ax2.spines[spine].set_visible(False)
    
    ax2.fill_between(np.arange(len(P_NewSigs)),P_NewSigs,color = 'b')
#    index,sigs = UpdateDI(NewSigs)
#    ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
#    ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
    ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    ax2.set_xlim(0,len(P_NewSigs))
    ax2.set_ylim(0,min(M_NewSigs.max(),P_NewSigs.max()))
    plt.savefig(P_out)
    plt.close(fig2)



def Plotting_SubArea_With_Sigs_Alleic_New2(M_out,P_out,M_Lib,P_Lib,M_Sigs,P_Sigs,
                                   chro,start,end,res,ls,lm,le):
    """
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    M_Matrix = M_Lib[chro]
    P_Matrix = P_Lib[chro]
    s = start // res
    e = end // res
    M_M = M_Matrix[s:e,s:e]
    P_P = P_Matrix[s:e,s:e]
    M_NewSigs = M_Sigs[s*40:e*40]
    P_NewSigs = P_Sigs[s*40:e*40]
    ls_ = ls // 1000 - start // 1000
    lm_ = lm // 1000 - start // 1000
    le_ = le // 1000 - start // 1000
    M_NewSigs[:ls_-40] = 0
    M_NewSigs[ls_+80:lm_-40] = 0
    M_NewSigs[lm_+80:le_-40] = 0
    M_NewSigs[le_+40:] = 0
    P_NewSigs[:ls_-40] = 0
    P_NewSigs[ls_+80:lm_-40] = 0
    P_NewSigs[lm_+80:le_-40] = 0
    P_NewSigs[le_+40:] = 0
    
    nonzero = M_M[np.nonzero(M_M)]
    #vmax = np.percentile(nonzero,90)
    size = (10,9)
    width = 0.618;Left = (1-width) / 2
    HB = 0.1;HH = width * size[0] / size[1]
    SB = HB + HH
    
    fig1 = plt.figure(figsize = size)
    ax = fig1.add_axes([Left,HB,width,HH])
    ax.imshow(M_M,cmap = my_cmap, aspect = 'auto',interpolation = 'none',
              extent = (0,M_M.shape[0],0,M_M.shape[0]),vmax = 37, origin = 'lower')
    
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    ax2 = fig1.add_axes([Left,SB,width,0.1])
    for spine in ['right','top','left']:
        ax2.spines[spine].set_visible(False)
    
    ax2.fill_between(np.arange(len(M_NewSigs)),M_NewSigs,color = 'b')
#    index,sigs = UpdateDI(NewSigs)
#    ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
#    ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
    ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    ax2.set_xlim(0,len(M_NewSigs))
    ax2.set_ylim(0,max(M_NewSigs.max(),P_NewSigs.max())*0.9)
    plt.savefig(M_out)
    plt.close(fig1)


    fig2 = plt.figure(figsize = size)
    ax = fig2.add_axes([Left,HB,width,HH])
    ax.imshow(P_P,cmap = my_cmap, aspect = 'auto',interpolation = 'none',
              extent = (0,P_P.shape[0],0,P_P.shape[0]),vmax = 37, origin = 'lower')
    
    ax.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    
    ax2 = fig2.add_axes([Left,SB,width,0.1])
    for spine in ['right','top','left']:
        ax2.spines[spine].set_visible(False)
    
    ax2.fill_between(np.arange(len(P_NewSigs)),P_NewSigs,color = 'b')
#    index,sigs = UpdateDI(NewSigs)
#    ax2.fill_between(index,sigs,where=sigs<=0,color = '#7093DB')
#    ax2.fill_between(index,sigs,where=sigs>=0,color = '#E47833')
    ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
    ax2.set_xlim(0,len(P_NewSigs))
    ax2.set_ylim(0,min(M_NewSigs.max(),P_NewSigs.max()) * 0.9)
    plt.savefig(P_out)
    plt.close(fig2)



def Plotting_with_Sigs(prefix,Lib,Sigs,res,interval,data_type = 'A'):
    """
    """
    #-------settings-----------
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    plot_sets = plotting_settings(ResHiC = res,interval = interval)
#    chrom = plot_sets['chroms']
    ResHiC = plot_sets['ResHiC']
    interval = plot_sets['interval']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]
    Step = int(res * interval / 1000)
    #------DataLoading---------
    #HiC_data = Load_HiC_Data(npz)
    
    #------Drawing---------
    chrom = Lib.keys()
    for i in chrom:
        inter  = Lib[i]
        pp = PdfPages(prefix + '-' + str(i) + '.pdf')
        startHiC = 0
        startSig = 0
        if data_type == 'A':
            Len = inter['bin1'].max() // interval
        else:
            Len = inter.shape[0] // interval
        for idx in range(Len):
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left,HB,width,HH])
            EndHiC = startHiC + interval
            EndSig = startSig + Step
            
            tmp_Sigs = Sigs[i][startSig:EndSig]
            Matrix = getmatrix(inter,startHiC,EndHiC,data_type)
            nonzero = Matrix[np.nonzero(Matrix)]
            if nonzero.size <= 100:
                plt.close(fig)
                startHiC = EndHiC
                startSig += Step
                continue
            
            vmax = np.percentile(nonzero,95)
            
            sc = ax.imshow(Matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, interval, 0, interval), vmax = vmax, origin = 'lower')
            
            ticks = list(np.linspace(0,interval,5).astype(int))
            pos = [((startHiC + t) * ResHiC) for t in ticks]
            labels = [properU(p) for p in pos]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            
            ax = fig.add_axes([Left + width + 0.02, HB, 0.01, HH])
            fig.colorbar(sc,cax = ax)
            
            ax2 = fig.add_axes([Left,HB+HH,width,0.1])
            for spine in ['right','top']:
                ax2.spines[spine].set_visible(False)
            
            ax2.fill_between(np.arange(len(tmp_Sigs)),tmp_Sigs,color = 'red')
            ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                            right = False, labelbottom = False, labeltop = False,
                            labelleft = False, labelright = False)
            ax2.set_xlim(0,(EndSig - startSig))
            ax2.set_ylim(0,tmp_Sigs.max())
            startHiC = EndHiC
            startSig = EndSig
            pp.savefig(fig)
            plt.close(fig)
        pp.close()   



def Load_Peaks(peak_fil):
    """
    """
    f = open(peak_fil,'r')
    Peaks = []
    for line in f:
        line = line.strip().split()
        tmp = (line[0].lstrip('chr'),line[1],line[2],line[6])
        Peaks.append(tmp)
    
    dtype = np.dtype({'names':['chr','start','end','value'],
                      'formats':['S4',np.int,np.int,np.float]})
    
    Peaks = np.array(Peaks,dtype = dtype)
    
    return Peaks



def Plotting_with_Peaks(prefix,Lib,peaks,res,interval,data_type = 'A'):
    """
    """
    #-------settings-----------
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    plot_sets = plotting_settings(ResHiC = res,interval = interval)
#    chrom = plot_sets['chroms']
    ResHiC = plot_sets['ResHiC']
    interval = plot_sets['interval']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]
    #------DataLoading---------
    #HiC_data = Load_HiC_Data(npz)
    
    #------Drawing---------
    #chrom = Lib.keys()
    chrom = ['1']
    for i in chrom:
        print i
        inter  = Lib[i]
        pp = PdfPages(prefix + '-' + str(i) + '.pdf')
        startHiC = 0
        if data_type == 'A':
            Len = inter['bin1'].max() // interval
        else:
            Len = inter.shape[0] // interval
        
        Chro_Peaks = peaks[peaks['chr'] == i]
        for idx in range(Len):
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left,HB,width,HH])
            EndHiC = startHiC + interval

            Matrix = getmatrix(inter,startHiC,EndHiC,data_type)
            nonzero = Matrix[np.nonzero(Matrix)]
            if nonzero.size <= 100:
                plt.close(fig)
                startHiC = EndHiC
                continue
            
            start = startHiC * res
            end = EndHiC * res
                        
            mask = (Chro_Peaks['start'] >= start) & (Chro_Peaks['end'] <= end)
            tmp_peak = Chro_Peaks[mask]
            peak_Sigs = {}
            for cp in tmp_peak:
                index = (cp['start'] + cp['end']) // 40000
                try:
                    peak_Sigs[index] += cp['value']
                except:
                    peak_Sigs[index] = 0
                    peak_Sigs[index] += cp['value']
            
            index = np.array(peak_Sigs.keys()) - startHiC
            values = peak_Sigs.values()
            
            vmax = np.percentile(nonzero,95)
            sc = ax.imshow(Matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, interval, 0, interval), vmax = vmax, origin = 'lower')
            
            ticks = list(np.linspace(0,interval,5).astype(int))
            pos = [((startHiC + t) * ResHiC) for t in ticks]
            labels = [properU(p) for p in pos]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            
            ax = fig.add_axes([Left + width + 0.02, HB, 0.01, HH])
            fig.colorbar(sc,cax = ax)
            
            ax2 = fig.add_axes([Left,HB+HH,width,0.1])
            for spine in ['right','top']:
                ax2.spines[spine].set_visible(False)
            ax2.bar(index,values,color = 'red')
            ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                            right = False, labelbottom = False, labeltop = False,
                            labelleft = False, labelright = False)

            startHiC = EndHiC
            pp.savefig(fig)
            plt.close(fig)
        pp.close()


def Plotting_whole_chromsome(prefix,Lib,res):
    """
    """
    Matplotlib_setting()
    my_cmap = Cmap_setting()
    plot_sets = plotting_settings(ResHiC = res)
    ResHiC = plot_sets['ResHiC']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]
    chrom = Lib.keys()
    pp = PdfPages(prefix+'.pdf')
    for i in chrom:
        Matrix = Lib[i]
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left,HB,width,HH])
        nonzero = Matrix[np.nonzero(Matrix)]
        
        vmax = np.percentile(nonzero,98)
        
        sc = ax.imshow(Matrix,cmap=my_cmap,aspect = 'auto',interpolation = 'none',
                       extent = (0,Matrix.shape[0],0,Matrix.shape[0]), vmax = vmax,
                        origin = 'lower')
        
        ticks = list(np.linspace(0,Matrix.shape[0],5).astype(int))
        pos = [(t * ResHiC) for t in ticks]
        labels = [properU(p) for p in pos]
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_title('chr'+str(i))
        ax = fig.add_axes([Left + width + 0.02, HB, 0.01, HH])
        fig.colorbar(sc,cax = ax)
        pp.savefig(fig)
        plt.close(fig)
    
    pp.close()
    

            
def Plotting_Loop(prefix,npz,cluster_loop,data_type = 'B'):
    """
        Plotting ...
        Parameters
        ----------
        prefix : str
            prefix of Output pdfs.
        
        npz : str
            HiC data.(npz)
        
        data_type : str
            A : 3 tuples Matrix for each chromosome in npz dict.
            B : whole 2D Matrix for each chromosome in npz dict.
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    plot_sets = plotting_settings()
    ResHiC = plot_sets['ResHiC']
    #interval = plot_sets['interval']
    interval = 100
    #chrom = plot_sets['chroms']
    #chrom = ['1']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]

    #------DataLoading---------
    #HiC_data = Load_HiC_Data(npz)
    
    #------Drawing---------
    #raw_loop = Load_loops(raw_loop)
    Cluster_loop = Load_Addtional_Loop(cluster_loop)
    Lib = np.load(npz)
    Lib = Lib['Matrix'][()]
    chrom = Lib.keys()
    for i in chrom:
        inter  = Lib[i]
        #raw_loops = raw_loop[raw_loop['chr'] == i]
        Cluster_loops = Cluster_loop[Cluster_loop['chr'] == i]
        pp = PdfPages(prefix + '-' + str(i) + '.pdf')
        startHiC = 0
        if data_type == 'A':
            Len = inter['bin1'].max() // interval
        else:
            Len = inter.shape[0] // interval
        for idx in range(Len):
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left,HB,width,HH])
            EndHiC = startHiC + interval
            
            Matrix = getmatrix(inter,startHiC,EndHiC,data_type)
            nonzero = Matrix[np.nonzero(Matrix)]
            if nonzero.size < 100:
                plt.close(fig)
                startHiC = EndHiC
                continue
            vmax = np.percentile(nonzero,96)
#            raw_mask = (raw_loops['end'] < EndHiC * ResHiC) & \
#                        (raw_loops['start'] >= startHiC * ResHiC)
            cluster_mask = (Cluster_loops['end'] < EndHiC * ResHiC) & \
                        (Cluster_loops['start'] >= startHiC * ResHiC)
                        
#            raw_extract_loop = raw_loops[raw_mask]
            cluster_extract_loop = Cluster_loops[cluster_mask]
            sc = ax.imshow(Matrix, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0, interval, 0, interval), vmax = vmax, origin = 'lower')
#            for l in raw_extract_loop:
#                x = l['start'] // ResHiC - startHiC
#                y = l['end'] // ResHiC - startHiC
#                ax.scatter(x+0.5,y+0.5,color = '',edgecolors = 'b', s = 10)

            for l in cluster_extract_loop:
                x = l['start'] // ResHiC - startHiC
                y = l['end'] // ResHiC - startHiC
                ax.scatter(y+0.5,x+0.5,color = '',edgecolors = 'b', s = 10)
                ax.text(y-2.5,x-3.5,'%s' % l['P_V'])
                ax.text(y-2.5,x-5.5,'%s' % l['Sig'])
                
            ticks = list(np.linspace(0,interval,5).astype(int))
            pos = [((startHiC + t) * ResHiC) for t in ticks]
            labels = [properU(p) for p in pos]
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            ax.set_xlim(0,interval)
            ax.set_ylim(0,interval)
            
            ax = fig.add_axes([Left + width + 0.02, HB, 0.01, HH])
            fig.colorbar(sc,cax = ax)
            
            startHiC = EndHiC
            pp.savefig(fig)
            plt.close(fig)
        pp.close()    





def Test_Plot(Matrix,out,Loop,Cluster_Loop,Gap):
    """
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    ResHiC = 40000
    interval = 200
    size = (12,14)
    
    Len = Matrix.shape[0] // interval
    
    pp = PdfPages(out+'.pdf')
    startHiC = 0
    raw_loops = Load_loops(Loop)
    Cluster_loops = Load_cluster_loops(Cluster_Loop)
#    Loop_Strength = []
    for idx in range(Len):
        fig, ax = plt.subplots(figsize = size)
        EndHiC = startHiC + interval
        
        M = Matrix[startHiC:EndHiC,startHiC:EndHiC]
        nonzero = M[np.nonzero(M)]
        raw_mask = (raw_loops['end'] < EndHiC * ResHiC) & \
                        (raw_loops['start'] >= startHiC * ResHiC)
        cluster_mask = (Cluster_loops['end'] < EndHiC * ResHiC) & \
                        (Cluster_loops['start'] >= startHiC * ResHiC)
                        
        raw_extract_loop = raw_loops[raw_mask]
        cluster_extract_loop = Cluster_loops[cluster_mask]
        if nonzero.size <= 100:
            plt.close(fig)
            startHiC = EndHiC
            continue

        mask = (Gap['start'] >= startHiC) & (Gap['end'] < EndHiC)
        sub_gap = Gap[mask]
        
        vmax = np.percentile(nonzero,95)
        sc = ax.imshow(M,cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0,interval,0,interval),vmax = vmax, origin = 'lower')
                       
        for l in raw_extract_loop:
            x = l['start'] // ResHiC - startHiC
            y = l['end'] // ResHiC - startHiC
            ax.scatter(x+0.5,y+0.5,color = '',edgecolors = 'b', s = 10)

        for l in cluster_extract_loop:
            x = l['start'] // ResHiC - startHiC
            y = l['end'] // ResHiC - startHiC
            ax.scatter(y+0.5,x+0.5,color = '',edgecolors = 'b', s = 10)
            ax.text(y+3.5,x-2.5,'%s' % l['wq'])
#            Loop_Strength.append(M[x][y])
        ticks = list(np.linspace(0,interval,5).astype(int))
        pos = [((startHiC + t) * ResHiC) for t in ticks]
        labels = [properU(p) for p in pos]
        
        for gap in sub_gap:
#            ax.plot([min(ticks),max(ticks)],[gap['start']-startHiC,gap['start']-startHiC],color='blue',
#                     lw = 0.5, ls = ':')
#            
#            ax.plot([min(ticks),max(ticks)],[gap['end']-startHiC,gap['end']-startHiC],color='green',
#                     lw = 0.5, ls = ':')
            start = gap['start'] - startHiC
            end = gap['end'] - startHiC
            if start == 0:
                ax.plot([0,end],[end,end], color = '#000000',
                        lw = 0.75)
                ax.plot([end,end],[0,end], color = '#000000',
                        lw = 0.75)
            elif end == interval:
                ax.plot([start,start],[start,interval],color = '#000000',
                        lw = 0.75)
                ax.plot([start,interval],[start,start],color = '#000000',
                        lw = 0.75)
            else:
                ax.plot([start,start],[start,end],
                        color = '#000000', lw = 0.75)
                ax.plot([start,end],[start,start],
                        color = '#000000', lw = 0.75)
                ax.plot([start,end],[end,end],
                        color = '#000000', lw = 0.75)
                ax.plot([end,end],[start,end],
                        color = '#000000', lw = 0.75)
                      
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_xlim(min(ticks),max(ticks))
        ax.set_ylim(min(ticks),max(ticks))
        
    
        fig.colorbar(sc)
        
        startHiC = EndHiC
        pp.savefig(fig)
        plt.close(fig)
    
    pp.close()
#    return Loop_Strength

def Whole_Chro_Plot(NPZ,out,cluster_loops,raw_loops):
    """
    """
    Matplotlib_setting()
    my_cmap  = Cmap_setting()
    plot_sets = plotting_settings(ResHiC = 20000)
    ResHiC = plot_sets['ResHiC']
    interval = plot_sets['interval']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]
    
#    Len = Matrix.shape[0] // interval
    
    #Loading Data
    Lib_M = np.load(NPZ)
    Lib_M = Lib_M['Matrix'][()]
    cluster_loops = Load_loops(cluster_loops)
    raw_loops = Load_loops(raw_loops)    
    for c in Lib_M.keys():
        print "chromosome %s" % c
        pp = PdfPages(out+'_'+c+'_'+'.pdf')
        Matrix = Lib_M[c]
        Len = Matrix.shape[0] // interval
        startHiC = 0        
        for idx in range(Len):
            fig = plt.figure(figsize = size)
            ax = fig.add_axes([Left,HB,width,HH])
            EndHiC = startHiC + interval
            
            M = Matrix[startHiC:EndHiC,startHiC:EndHiC]
            nonzero = M[np.nonzero(M)]
            
            if nonzero.size <= 100:
                plt.close(fig)
                startHiC = EndHiC
                continue
            raw_mask = (raw_loops['end'] < EndHiC * ResHiC) & \
                        (raw_loops['start'] >= startHiC * ResHiC)
            raw_extract_loop = raw_loops[raw_mask]
        
            cluster_mask = (cluster_loops['end'] < EndHiC * ResHiC) & \
                           (cluster_loops['start'] >= startHiC * ResHiC)
            cluster_extract_loop = cluster_loops[cluster_mask]
            
            vmax = np.percentile(nonzero,95)
            sc = ax.imshow(M,cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                           extent = (0,interval,0,interval),vmax = vmax, origin = 'lower')
            ticks = list(np.linspace(0,interval,5).astype(int))
            pos = [((startHiC + t) * ResHiC) for t in ticks]
            labels = [properU(p) for p in pos]            
            
            for l in raw_extract_loop:
                x = l['start'] // ResHiC - startHiC
                y = l['end'] // ResHiC - startHiC
                ax.scatter(x+0.5,y+0.5,color = '',edgecolors = 'b', s = 10)
                
            for l in cluster_extract_loop:
                x = l['start'] // ResHiC - startHiC
                y = l['end'] // ResHiC - startHiC
                ax.scatter(y+0.5,x+0.5,color = '', edgecolors = 'b', s = 10)
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            ax.set_xlim(min(ticks),max(ticks))
            ax.set_ylim(min(ticks),max(ticks))
            
            ax = fig.add_axes([Left+width+0.02,HB,0.01,HH])
            fig.colorbar(sc,cax = ax)
            
            startHiC = EndHiC
            pp.savefig(fig)
            plt.close(fig)
        
        pp.close()
            
            
            
            
def Allel_AP_HeatMap(Matrix,out,cluster_loop):
    """
    """
    Matplotlib_setting()
    #my_cmap  = Cmap_setting(start_Color = '#0000FF',end_Color = '#CD0000')
    my_cmap = LinearSegmentedColormap.from_list('interactions',['#0000FF','#FFFFFF','#CD0000'])
    plot_sets = plotting_settings(ResHiC = 20000)
    ResHiC = plot_sets['ResHiC']
    interval = plot_sets['interval']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]
    
    
    
    Len = Matrix.shape[0] // interval
    
    pp = PdfPages(out+'.pdf')
    startHiC = 0    
    for idx in range(Len):
        fig = plt.figure(figsize = size)
        ax = fig.add_axes([Left,HB,width,HH])
        EndHiC = startHiC + interval
        
        M = Matrix[startHiC:EndHiC,startHiC:EndHiC]
        nonzero = M[np.nonzero(M)]
        
        if nonzero.size <= 100:
            plt.close(fig)
            startHiC = EndHiC
            continue
        
        
        sc = ax.imshow(M,cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0,interval,0,interval),vmin = -4,vmax = 4, origin = 'lower')
        
        ticks = list(np.linspace(0,interval,5).astype(int))
        pos = [((startHiC + t) * ResHiC) for t in ticks]
        labels = [properU(p) for p in pos]
        
            
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_xlim(min(ticks),max(ticks))
        ax.set_ylim(min(ticks),max(ticks))
        ax = fig.add_axes([Left + width + 0.02, HB, 0.01, HH])
        fig.colorbar(sc,cax = ax)
        
        startHiC = EndHiC
        pp.savefig(fig)
        plt.close(fig)
    
    pp.close()
    





def Allelic_HeatMap_Loops(M_npz,P_npz,loop,out):
    """
    """
    #Loading Data
    M_Lib = np.load(M_npz)
    P_Lib = np.load(P_npz)
    loop_type = np.dtype({'names':['chr','start','end','M_C','P_C'],
                          'formats':['S4',np.int,np.int,np.int,np.int]})
    loops = np.loadtxt(loop,dtype = loop_type, usecols = (0,1,2,3,4))
    
    # PDF setting   
    pp_M = PdfPages(out+'_Maternal.pdf')
    pp_P = PdfPages(out+'_Paternal.pdf')
    Matplotlib_setting()
    my_cmap  = Cmap_setting(start_Color = '#FFFFFF',end_Color = '#CD0000')
    #my_cmap = LinearSegmentedColormap.from_list('interactions',['#0000FF','#FFFFFF','#CD0000'])
    plot_sets = plotting_settings(ResHiC = 20000)
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]
    
    count = 0
    for l in loops:
        count += 1
        print count
        chro = l['chr']
        s = (l['start'] - 2000000) // 20000
        e = (l['end'] + 2000000) // 20000
        
        M_M = M_Lib[chro][s:e,s:e]
        P_M = P_Lib[chro][s:e,s:e]
        
        M_nonzero = M_M[np.nonzero(M_M)]
        P_nonzero = P_M[np.nonzero(P_M)]
        
        if M_nonzero.size <= 100 or P_nonzero.size <= 100:
            plt.close(fig_M)
            plt.close(fig_P)
            continue       
        vmax = np.percentile(M_nonzero,95)
        ticks = [0,e-s]
        labels = [properU((p+s)*20000) for p in ticks]
        shape = M_M.shape[0]
        
        #Maternal
        fig_M = plt.figure(figsize = size)
        ax = fig_M.add_axes([Left,HB,width,HH])
        
        M_HM = ax.imshow(M_M.T, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                        extent = (0,shape,0,shape),vmax = vmax, origin = 'lower')
        
        x = l['start'] // 20000 - s
        y = l['end'] // 20000 - s
        ax.scatter(x+0.5,y+0.5,color = '', edgecolor = 'b', s = 10)
        ax.scatter(y+0.5,x+0.5,color = '', edgecolor = 'b', s = 10)
        ax.text(y-1.5,x-7.5,'%s' % l['M_C'],size = 15)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_xlim(min(ticks),max(ticks))
        ax.set_ylim(min(ticks),max(ticks))
        title = l[0]+':  '+properU(l[1])+'-----'+properU(l[2])
        ax.set_title(title)
        ax = fig_M.add_axes([Left + width + 0.02, HB, 0.01, HH])
        fig_M.colorbar(M_HM,cax = ax)
        
        pp_M.savefig(fig_M)
        plt.close(fig_M)
        
        #Paternal
        fig_P = plt.figure(figsize = size)
        ax = fig_P.add_axes([Left,HB,width,HH])
        
        P_HM = ax.imshow(P_M.T, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                         extent = (0,shape,0,shape),vmax = vmax, origin = 'lower')
        
        x = l['start'] // 20000 - s
        y = l['end'] // 20000 - s
        ax.scatter(x+0.5,y+0.5,color = '', edgecolor = 'b', s = 10)
        ax.scatter(y+0.5,x+0.5,color = '', edgecolor = 'b', s = 10)
        
        ax.text(y-1.5,x-7.5,'%s' % l['P_C'],size = 15)
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels)
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)
        ax.set_xlim(min(ticks),max(ticks))
        ax.set_ylim(min(ticks),max(ticks))
        title = l[0]+':  '+properU(l[1])+'-----'+properU(l[2])
        ax.set_title(title)
        ax = fig_P.add_axes([Left + width + 0.02, HB, 0.01, HH])
        fig_P.colorbar(P_HM,cax = ax)
        
        pp_P.savefig(fig_P)
        plt.close(fig_P)
    
    pp_M.close()
    pp_P.close()



def Load_Loop_Select(loop_fil):
    """
    """
    loop_type = np.dtype({'names':['chr','start','end'],
                          'formats':['S4',np.int,np.int]})
    
    loops = np.loadtxt(loop_fil, dtype = loop_type, skiprows = 1, usecols = [0,1,2])
    
    return loops



def Loop_stay(Matrix,loop_x,loop_y):
    """
    """
    #for loci (x,y)
    Mark_loci_1 = False
    
    Left_Inter_1 = Matrix[loop_x,loop_y-3:loop_y]
    Right_Inter_1 = Matrix[loop_x,loop_y+1:loop_y+4]
    Top_Inter_1 = Matrix[loop_x-3:loop_x,loop_y]
    Bottom_Inter_1 = Matrix[loop_x+1:loop_x+4,loop_y]
    
    Left_Zero_1 = (Left_Inter_1==0).sum()
    Right_Zero_1 = (Right_Inter_1==0).sum()
    Top_Zero_1 = (Top_Inter_1==0).sum()
    Bottom_Zero_1 = (Bottom_Inter_1==0).sum()
    
    
    if Left_Zero_1 > 0 and Right_Zero_1 > 0  and (Left_Zero_1 + Right_Zero_1) > 2:
        pass
    elif Top_Zero_1 > 0 and Bottom_Zero_1 > 0 and (Top_Zero_1 + Bottom_Zero_1) > 2:
        pass
    else:
        Mark_loci_1 = True
    
    #for loci (y,x)
    Mark_loci_2 = False
    
    Left_Inter_2 = Matrix[loop_y,loop_x-3:loop_x]
    Right_Inter_2 = Matrix[loop_y,loop_x+1:loop_x+4]
    Top_Inter_2 = Matrix[loop_y-3:loop_y,loop_x]
    Bottom_Inter_2 = Matrix[loop_y+1:loop_y+4,loop_x]
    
    Left_Zero_2 = (Left_Inter_2==0).sum()
    Right_Zero_2 = (Right_Inter_2==0).sum()
    Top_Zero_2 = (Top_Inter_2==0).sum()
    Bottom_Zero_2 = (Bottom_Inter_2==0).sum()
    
    if Left_Zero_2 > 0 and Right_Zero_2 > 0  and (Left_Zero_2 + Right_Zero_2) > 2:
        pass
    elif Top_Zero_2 > 0 and Bottom_Zero_2 > 0 and (Top_Zero_2 + Bottom_Zero_2) > 2:
        pass
    else:
        Mark_loci_2 = True
    
    
    #for loop
    if Mark_loci_1 or Mark_loci_2:
        return True
    else:
        return False


def Loop_Select(NPZ_M,NPZ_P,loop_fil,out):
    """
    """
    Lib_M = np.load(NPZ_M)
    Lib_P = np.load(NPZ_P)
    f = open(loop_fil)
    c = ''
    with open(out,'w') as o:
        for line in islice(f,1,None):
            line = line.split()
            if line[0] != c:
                c = line[0]
                M_M = Lib_M[c]
                P_P = Lib_P[c]
            start = int(line[1]) // 20000
            end = int(line[2]) // 20000
            Loop_stay_M = Loop_stay(M_M,start,end)
            Loop_stay_P = Loop_stay(P_P,start,end)
            
            if Loop_stay_M or Loop_stay_P:
                o.writelines('\t'.join(line)+'\n')
    f.close()
    




  
def Loop_Select_plot(Matrix_M,Matrix_P,cluster_loop,out):
    """
    """
    Matplotlib_setting()
    my_cmap = Cmap_setting()
    plot_sets = plotting_settings(ResHiC = 20000)
    ResHiC = plot_sets['ResHiC']
    interval = plot_sets['interval']
    size = plot_sets['figsize']
    width = plot_sets['width']
    Left = (1 - width) / 2
    HB = plot_sets['HB']
    HH = width * size[0] / size[1]
    
    Len = Matrix_M.shape[0] // interval   
    
    
    cluster_loops = Load_Loop_Select(cluster_loop)

    cluster_loops = cluster_loops[cluster_loops['chr'] == '12']
    
    pp1 = PdfPages(out+'_Maternal.pdf')
    pp2 = PdfPages(out+'_Paternal.pdf')
    startHiC = 0 
    for idx in range(Len):
        fig1 = plt.figure(figsize = size)
        fig2 = plt.figure(figsize = size)
        ax1 = fig1.add_axes([Left,HB,width,HH])
        ax2 = fig2.add_axes([Left,HB,width,HH])
        EndHiC = startHiC + interval
        
        M_M = Matrix_M[startHiC:EndHiC,startHiC:EndHiC]
        P_P = Matrix_P[startHiC:EndHiC,startHiC:EndHiC]
        
        nonzero_M = M_M[np.nonzero(M_M)]
        nonzero_P = P_P[np.nonzero(P_P)]
        
        if nonzero_M.size <= 100 or nonzero_P.size <= 100:
            plt.close(fig1)
            plt.close(fig2)
            startHiC = EndHiC
            continue
        
        cluster_mask = (cluster_loops['end'] < EndHiC * ResHiC) & \
                       (cluster_loops['start'] >= startHiC * ResHiC)
                    
        cluster_extract_loop = cluster_loops[cluster_mask]
        
        vmax = min(np.percentile(nonzero_M,95),np.percentile(nonzero_P,95))
        sc1 = ax1.imshow(M_M,cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                       extent = (0,interval,0,interval),vmax = vmax, origin = 'lower')
        
        sc2 = ax2.imshow(P_P,cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                         extent = (0,interval,0,interval),vmax = vmax, origin = 'lower')
                         
        ticks = list(np.linspace(0,interval,5).astype(int))
        pos = [((startHiC + t) * ResHiC) for t in ticks]
        labels = [properU(p) for p in pos]
        
        for l in cluster_extract_loop:
            x = l['start'] // ResHiC - startHiC
            y = l['end'] // ResHiC - startHiC
            
            if (M_M[x,y] + M_M[y,x]) == 0 or (P_P[x,y] + P_P[y,x]) == 0:
                ax1.scatter(x+0.5,y+0.5,color = '',edgecolors = 'b', s = 10)
                ax1.scatter(y+0.5,x+0.5,color = '',edgecolors = 'b', s = 10)
                ax2.scatter(x+0.5,y+0.5,color = '',edgecolors = 'b', s = 10)
                ax2.scatter(y+0.5,x+0.5,color = '',edgecolors = 'b', s = 10)
            else:
                loop_stay_M = Loop_stay(M_M,x,y)
                loop_stay_P = Loop_stay(P_P,x,y)
                if loop_stay_M or loop_stay_P:
                    ax1.scatter(x+0.5,y+0.5,color = '',edgecolors = 'green', s = 10)
                    ax1.scatter(y+0.5,x+0.5,color = '',edgecolors = 'green', s = 10)
                    ax2.scatter(x+0.5,y+0.5,color = '',edgecolors = 'green', s = 10)
                    ax2.scatter(y+0.5,x+0.5,color = '',edgecolors = 'green', s = 10)
                else:
                    ax1.scatter(x+0.5,y+0.5,color = '',edgecolors = 'b', s = 10)
                    ax1.scatter(y+0.5,x+0.5,color = '',edgecolors = 'b', s = 10)
                    ax2.scatter(x+0.5,y+0.5,color = '',edgecolors = 'b', s = 10)
                    ax2.scatter(y+0.5,x+0.5,color = '',edgecolors = 'b', s = 10)
            
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labels)
        ax1.set_yticks(ticks)
        ax1.set_yticklabels(labels)
        ax1.set_xlim(min(ticks),max(ticks))
        ax1.set_ylim(min(ticks),max(ticks))
        
        ax2.set_xticks(ticks)
        ax2.set_xticklabels(labels)
        ax2.set_yticks(ticks)
        ax2.set_yticklabels(labels)
        ax2.set_xlim(min(ticks),max(ticks))
        ax2.set_ylim(min(ticks),max(ticks))
        
        ax = fig1.add_axes([Left + width + 0.02, HB, 0.01, HH])
        fig1.colorbar(sc1,cax = ax)
        ax = fig2.add_axes([Left + width + 0.02, HB, 0.01, HH])
        fig2.colorbar(sc2,cax = ax)
        startHiC = EndHiC
        pp1.savefig(fig1)
        pp2.savefig(fig2)
        plt.close(fig1)
        plt.close(fig2)
    
    pp1.close()          
    pp2.close()

def Load_Sigs(fil):
    """
    """
    Sigs = {}
    f = open(fil,'r')
    for line in f:
        line = line.strip().split()
        if line[0] not in Sigs.keys():
            Sigs[line[0]] = []
            Sigs[line[0]].append((line[1],line[2],line[3]))
        else:
            Sigs[line[0]].append((line[1],line[2],line[3]))
    Sigs_type = np.dtype({'names':['start','end','value'],
                          'formats':[np.int,np.int,np.int]})
    
    for c, v in Sigs.items():
        v = np.array(v,dtype = Sigs_type)
        Sigs[c] = v
    
    return Sigs

def Get_Sigs_Values(loop,M_RNA_S,P_RNA_S,M_ATAC_S,P_ATAC_S):
    """
    """
    M_RNA_Values = []
    P_RNA_Values = []
    M_ATAC_Values = []
    P_ATAC_Values = []
    chro = 'chr'+loop['chr']
    start = loop['start'] - 200000
    end = loop['end'] + 200000
    index = np.arange(start,end+1000,1000)
    for i in range(len(index)-1):
        mask = (M_RNA_S['chr'] == chro) & \
                (M_RNA_S['start']<=index[i]) & (M_RNA_S['end']>=index[i+1])
        tmp = M_RNA_S[mask]
        if tmp.size == 0:
            M_RNA_Values.append(0)
        else:
            M_RNA_Values.append(tmp['value'])
        
        mask = (P_RNA_S['chr'] == chro) & \
                (P_RNA_S['start']<=index[i]) & (P_RNA_S['end']>=index[i+1])
        tmp = P_RNA_S[mask]
        if tmp.size == 0:
            P_RNA_Values.append(0)
        else:
            P_RNA_Values.append(tmp['value'])
        
        mask = (M_ATAC_S['chr'] == chro) & \
                (M_ATAC_S['start']<=index[i]) & (M_ATAC_S['end']>=index[i+1])
        tmp = M_ATAC_S[mask]
        if tmp.size == 0:
            M_ATAC_Values.append(0)
        else:
            M_ATAC_Values.append(tmp['value'])
        
        mask = (P_ATAC_S['chr'] == chro) & \
                (P_ATAC_S['start']<=index[i]) & (P_ATAC_S['end']>=index[i+1])
        tmp = P_ATAC_S[mask]
        if tmp.size == 0:
            P_ATAC_Values.append(0)
        else:
            P_ATAC_Values.append(tmp['value'])

    M_RNA_Values = np.array(M_RNA_Values)
    M_RNA_Values = M_RNA_Values.reshape(M_RNA_Values.shape[0],)
    P_RNA_Values = np.array(P_RNA_Values)
    P_RNA_Values = P_RNA_Values.reshape(P_RNA_Values.shape[0],)
    M_ATAC_Values = np.array(M_ATAC_Values)
    M_ATAC_Values = M_ATAC_Values.reshape(M_ATAC_Values.shape[0],)
    P_ATAC_Values = np.array(P_ATAC_Values)
    P_ATAC_Values = P_ATAC_Values.reshape(P_ATAC_Values.shape[0],)
    
    return M_RNA_Values,P_RNA_Values,M_ATAC_Values,P_ATAC_Values 
    
def Allelic_Gene_Loops(M_NPZ,P_NPZ,M_ATAC_fil,P_ATAC_fil,
                       M_RNA_fil,P_RNA_fil,RNA_fil,Loop_fil,out):
    """
    """
    #Data_Loading
    Bedgraph_type = np.dtype({'names':['chr','start','end','value'],
                              'formats':['S6',np.int,np.int,np.int]})
    
    Loop_type = np.dtype({'names':['chr','start','end','CCS_M','CCS_P','fESC_M','fESC_P','NT5_M','NT5_P','NT6_M','NT6_P'],
                         'formats':['S4',np.int,np.int,np.int,np.int,np.int,np.int,np.int,np.int,np.int,np.int]})
    
    M_RNA_S = np.loadtxt(M_RNA_fil,dtype = Bedgraph_type)
    P_RNA_S = np.loadtxt(P_RNA_fil,dtype = Bedgraph_type)
    
    M_ATAC_S = np.loadtxt(M_ATAC_fil,dtype = Bedgraph_type)
    P_ATAC_S = np.loadtxt(P_ATAC_fil,dtype = Bedgraph_type)
    
    M_Lib = np.load(M_NPZ)
    P_Lib = np.load(P_NPZ)
    
    M_RNA, P_RNA, Non_RNA = Load_Speci_RNA(RNA_fil)
    
    Loop = np.loadtxt(Loop_fil,dtype = Loop_type)
    
    #About fig
    hexcolors = Dark2_8.hex_colors
    my_cmap = LinearSegmentedColormap.from_list('interaction',
                                            ['#FFFFFF','#CD0000'])
    size = (12,14)
    width = 0.618;Left = (1 - width) / 2
    HB = 0.1; HH = width * size[0] / size[1]
    SB = HB + HH
    ST = 0.82
    SH = (ST - SB) / 3
    pp1 = PdfPages(out+'_Maternal.pdf')
    pp2 = PdfPages(out+'_Paternal.pdf')
    loop_num = 1
    for loop in Loop:
        print loop_num
        loop_num += 1
        chro = loop['chr']
        s_ = (loop['start'] - 1000000) // 20000
        e_ = (loop['end'] + 1000000) // 20000
        
        M_M = M_Lib[chro][s_:e_,s_:e_]
        P_P = P_Lib[chro][s_:e_,s_:e_]
        
        M_nonzero = M_M[np.nonzero(M_M)]
        P_nonzero = P_P[np.nonzero(P_P)]
        
        M_RNA_Sig,P_RNA_Sig,M_ATAC_Sig,P_ATAC_Sig = Get_Sigs_Values(loop,M_RNA_S,P_RNA_S,
                                                                    M_ATAC_S,P_ATAC_S)
        
        if M_nonzero.size <= 50 and P_nonzero <= 50:
            plt.close(fig1)
            plt.close(fig2)
            continue
        
        vmax = max(np.percentile(M_nonzero,95),np.percentile(P_nonzero,95))
        ticks = [0,e_-s_]
        labels = [properU((p+s_)*20000) for p in ticks]
        shape = M_M.shape[0]
        
        #=============Maternal===============
        #=============HeatMap
        fig1 = plt.figure(figsize = size)
        ax1 = fig1.add_axes([Left,HB,width,HH])
        
        M_H = ax1.imshow(M_M.T, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                        extent = (0,shape,0,shape), vmax = vmax, origin = 'lower')
        
        x = loop['start'] // 20000 - s_
        y = loop['end'] // 20000 - s_
        ax1.scatter(x+0.5,y+0.5,color = '', edgecolor = 'b', s = 10)
        ax1.scatter(y+0.5,x+0.5,color = '', edgecolor = 'b', s = 10)
        ax1.text(y-1.5,x-7.5,'%s' % loop['M'],size = 15)
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labels)
        ax1.set_yticks(ticks)
        ax1.set_yticklabels(labels)
        ax1.set_xlim(min(ticks),max(ticks))
        ax1.set_ylim(min(ticks),max(ticks))
        ax1.set_xlabel('chr%s' % chro,size = 15)
        caxis_H(ax1)
        
        ax = fig1.add_axes([Left + width + 0.02, HB, 0.01, HH])
        fig1.colorbar(M_H, cax = ax)
        
        #===============Sigs
        #RNA
        cb = SB
        count = 0
        ax = fig1.add_axes([Left,cb,width,SH])
        ax.fill_between(np.arange(M_RNA_Sig.size),M_RNA_Sig, color = hexcolors[count],
                        edgecolor = 'none')
        ax.set_ylabel('M_RNA', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        ax.set_xlim(0,M_RNA_Sig.size)
        ax.set_ylim(0,M_RNA_Sig.max()*1.1)
        caxis_S(ax, hexcolors[count])
        
        cb += SH
        count += 1
        ax = fig1.add_axes([Left,cb,width,SH])
        ax.fill_between(np.arange(P_RNA_Sig.size),P_RNA_Sig, color = hexcolors[count],
                        edgecolor = 'none')
        ax.set_ylabel('P_RNA', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        ax.set_xlim(0,P_RNA_Sig.size)
        ax.set_ylim(0,M_RNA_Sig.max()*1.1)
        caxis_S(ax, hexcolors[count])
        #ATAC
        cb += SH
        count += 1
        ax = fig1.add_axes([Left,cb,width,SH])
        ax.fill_between(np.arange(M_ATAC_Sig.size),M_ATAC_Sig, color = hexcolors[count],
                        edgecolor = 'none')
        ax.set_ylabel('M_ATAC', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        ax.set_xlim(0,M_ATAC_Sig.size)
        ax.set_ylim(0,M_ATAC_Sig.max()*1.1)
        caxis_S(ax, hexcolors[count])
        
        cb += SH
        count += 1
        ax = fig1.add_axes([Left,cb,width,SH])
        ax.fill_between(np.arange(P_ATAC_Sig.size),M_ATAC_Sig, color = hexcolors[count],
                        edgecolor = 'none')
        ax.set_ylabel('P_ATAC', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        ax.set_xlim(0,P_ATAC_Sig.size)
        ax.set_ylim(0,M_ATAC_Sig.max()*1.1)
        caxis_S(ax, hexcolors[count])
        
        #Gene_Name
        cb += SH
        count += 1
        start = loop['start']
        end = loop['end']
        cis_dis = 100000
        trans_dis = 100000
        RNA_mask_M = (M_RNA['chr'] == chro) & \
                (((M_RNA['start'] >= start - trans_dis) & (M_RNA['end'] <= start + cis_dis)) | \
                ((M_RNA['start'] >= end - cis_dis) & (M_RNA['end'] <= end + trans_dis)))

        RNA_mask_P = (P_RNA['chr'] == chro) & \
                (((P_RNA['start'] >= start - trans_dis) & (P_RNA['end'] <= start + cis_dis)) | \
                ((P_RNA['start'] >= end - cis_dis) & (P_RNA['end'] <= end + trans_dis)))
        
        RNA_tmp_M = M_RNA[RNA_mask_M]
        RNA_tmp_P = P_RNA[RNA_mask_P]
        m_max = len(RNA_tmp_M)
        p_max = len(RNA_tmp_P)
        ax = fig1.add_axes([Left, cb, width, SH])
        ax.set_ylim(-p_max-2,m_max + 2)
        ax.set_xlim(min(ticks),max(ticks))
        
        m_N = 0 #Matenal Gene Count
        p_N = 0 #Patenal Gene Count
        for gene in RNA_tmp_M:
            m_N += 1
            x = (gene['start'] + gene['end']) // 2 // 20000 - s_ 
            ax.scatter(x,m_N,c = 'red',s = 8)
            ax.text(x+ 2,m_N,'%s' % gene['Name'],size = 8)
        for gene in RNA_tmp_P:
            p_N -= 1
            x = (gene['start'] + gene['end']) // 2 // 20000 - s_ 
            ax.scatter(x,p_N,c = 'blue',s = 8)
            ax.text(x+ 2,p_N,'%s' % gene['Name'],size = 8)
        
        ax.set_ylabel('Genes', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        caxis_S(ax,hexcolors[count])
        pp1.savefig(fig1)
        plt.close(fig1)
        
        #==============Paternal
        fig1 = plt.figure(figsize = size)
        ax1 = fig1.add_axes([Left,HB,width,HH])
        
        M_H = ax1.imshow(P_P.T, cmap = my_cmap, aspect = 'auto', interpolation = 'none',
                        extent = (0,shape,0,shape), vmax = vmax, origin = 'lower')
        
        s_ = (loop['start'] - 1000000) // 20000
        e_ = (loop['end'] + 1000000) // 20000
        
        x = loop['start'] // 20000 - s_
        y = loop['end'] // 20000 - s_
        ax1.scatter(x+0.5,y+0.5,color = '', edgecolor = 'b', s = 10)
        ax1.scatter(y+0.5,x+0.5,color = '', edgecolor = 'b', s = 10)
        ax1.text(y-1.5,x-7.5,'%s' % loop['P'],size = 15)
        ax1.set_xticks(ticks)
        ax1.set_xticklabels(labels)
        ax1.set_yticks(ticks)
        ax1.set_yticklabels(labels)
        ax1.set_xlim(min(ticks),max(ticks))
        ax1.set_ylim(min(ticks),max(ticks))
        ax1.set_xlabel('chr%s' % chro)
        caxis_H(ax1)
        
        ax = fig1.add_axes([Left + width + 0.02, HB, 0.01, HH])
        fig1.colorbar(M_H, cax = ax)
        
        #===============Sigs
        #RNA
        cb = SB
        count = 0
        ax = fig1.add_axes([Left,cb,width,SH])
        ax.fill_between(np.arange(M_RNA_Sig.size),M_RNA_Sig, color = hexcolors[count],
                        edgecolor = 'none')
        ax.set_ylabel('M_RNA', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        ax.set_xlim(0,M_RNA_Sig.size)
        ax.set_ylim(0,M_RNA_Sig.max()*1.1)
        caxis_S(ax, hexcolors[count])
        
        cb += SH
        count += 1
        ax = fig1.add_axes([Left,cb,width,SH])
        ax.fill_between(np.arange(P_RNA_Sig.size),P_RNA_Sig, color = hexcolors[count],
                        edgecolor = 'none')
        ax.set_ylabel('P_RNA', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        ax.set_xlim(0,P_RNA_Sig.size)
        ax.set_ylim(0,M_RNA_Sig.max()*1.1)
        caxis_S(ax, hexcolors[count])
        #ATAC
        cb += SH
        count += 1
        ax = fig1.add_axes([Left,cb,width,SH])
        ax.fill_between(np.arange(M_ATAC_Sig.size),M_ATAC_Sig, color = hexcolors[count],
                        edgecolor = 'none')
        ax.set_ylabel('M_ATAC', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        ax.set_xlim(0,M_ATAC_Sig.size)
        ax.set_ylim(0,M_ATAC_Sig.max()*1.1)
        caxis_S(ax, hexcolors[count])
        
        cb += SH
        count += 1
        ax = fig1.add_axes([Left,cb,width,SH])
        ax.fill_between(np.arange(P_ATAC_Sig.size),M_ATAC_Sig, color = hexcolors[count],
                        edgecolor = 'none')
        ax.set_ylabel('P_ATAC', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        ax.set_xlim(0,P_ATAC_Sig.size)
        ax.set_ylim(0,M_ATAC_Sig.max()*1.1)
        caxis_S(ax, hexcolors[count])
        
        #Gene_Name
        cb += SH
        count += 1
        start = loop['start']
        end = loop['end']
        cis_dis = 100000
        trans_dis = 100000
        RNA_mask_M = (M_RNA['chr'] == chro) & \
                (((M_RNA['start'] >= start - trans_dis) & (M_RNA['end'] <= start + cis_dis)) | \
                ((M_RNA['start'] >= end - cis_dis) & (M_RNA['end'] <= end + trans_dis)))

        RNA_mask_P = (P_RNA['chr'] == chro) & \
                (((P_RNA['start'] >= start - trans_dis) & (P_RNA['end'] <= start + cis_dis)) | \
                ((P_RNA['start'] >= end - cis_dis) & (P_RNA['end'] <= end + trans_dis)))
        
        RNA_tmp_M = M_RNA[RNA_mask_M]
        RNA_tmp_P = P_RNA[RNA_mask_P]
        m_max = len(RNA_tmp_M)
        p_max = len(RNA_tmp_P)
        ax = fig1.add_axes([Left, cb, width, SH])
        ax.set_ylim(-p_max-2,m_max + 2)
        ax.set_xlim(min(ticks),max(ticks))
        
        m_N = 0 #Matenal Gene Count
        p_N = 0 #Patenal Gene Count
        for gene in RNA_tmp_M:
            m_N += 1
            x = (gene['start'] + gene['end']) // 2 // 20000 - s_
            ax.scatter(x,m_N,c = 'red',s = 8)
            ax.text(x+ 2,m_N,'%s' % gene['Name'],size = 8)
        for gene in RNA_tmp_P:
            p_N -= 1
            x = (gene['start'] + gene['end']) // 2 // 20000 - s_
            ax.scatter(x,p_N,c = 'blue',s = 8)
            ax.text(x+ 2,p_N,'%s' % gene['Name'],size = 8)
        
        ax.set_ylabel('Genes', labelpad = 50, rotation = 'horizontal',
                              style = 'italic', size = 10)
        caxis_S(ax,hexcolors[count])
        pp2.savefig(fig1)
        plt.close(fig1)
    
    pp1.close()
    pp2.close()    


def Sax_Set(ax, color):
    """
    """
    for spine in ['right','top']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis = 'both', bottom = False, top = False, left = False,
                   right = False, labelbottom = False, labeltop = False,
                   labelleft = False, labelright = False)
    ax.spines['left'].set_lw(1.5)
    ax.spines['left'].set_color(color)
    ax.spines['left'].set_alpha(0.9)
    ax.spines['left'].set_linestyle('dotted')
    
    ax.spines['bottom'].set_lw(0.5)
    ax.spines['bottom'].set_color('#b3b3b3')
    ax.spines['bottom'].set_alpha(0.9)

    

    
def Plot_with_All_Sigs(out,res,Loops,
                       Lib,CTCF_Sigs,RNA_Sigs,ATAC_Sigs,
                       M_Lib,CTCF_Sigs_M,RNA_Sigs_M,ATAC_Sigs_M,
                       P_Lib,CTCF_Sigs_P,RNA_Sigs_P,ATAC_Sigs_P):
    """
    """
    Matplotlib_setting()
    my_cmap = Cmap_setting()
    ResHiC = res
    interval = 100
    Step = int(ResHiC * interval) // 1000
    size = (12,14)
    width = 0.618;Left = (1 - width) / 2
    HB = 0.1;HH = width * size[0] / size[1]
    SB = HB + HH
    ST = 0.86
    SH = (ST - SB) / 3
    for i in Lib.keys():
        pp_R = PdfPages(out+'Raw-'+str(i)+'.pdf')
        pp_M = PdfPages(out+'Maternal-'+str(i)+'.pdf')
        pp_P = PdfPages(out+'Paternal-'+str(i)+'.pdf')
        inter_R = Lib[i]
        inter_M = M_Lib[i]
        inter_P = P_Lib[i]
        startHiC = 0
        startSig = 0
        for idx in range(inter_R.shape[0] // interval):
            
            fig_R = plt.figure(figsize = size)
            fig_M = plt.figure(figsize = size)
            fig_P = plt.figure(figsize = size)
            
            axR1 = fig_R.add_axes([Left,HB,width,HH])
            axM1 = fig_M.add_axes([Left,HB,width,HH])
            axP1 = fig_P.add_axes([Left,HB,width,HH])
                        
            EndHiC = startHiC + interval
            EndSig = startSig + Step
            
            Matrix_R = getmatrix(inter_R,startHiC,EndHiC,data_type='B')
            Matrix_M = getmatrix(inter_M,startHiC,EndHiC,data_type='B')
            Matrix_P = getmatrix(inter_P,startHiC,EndHiC,data_type='B')
            
            nonzero_R = Matrix_R[np.nonzero(Matrix_R)]
            nonzero_M = Matrix_M[np.nonzero(Matrix_M)]
            nonzero_P = Matrix_P[np.nonzero(Matrix_P)]
            
            mask = (Loops['chr'] == i) & (Loops['start'] > startHiC * res) & (Loops['end'] < EndHiC * res)
            tmp_loops = Loops[mask]
            if tmp_loops.size == 0:
                plt.close(fig_R)
                plt.close(fig_M)
                plt.close(fig_P)
                startHiC = EndHiC
                startSig += Step
                continue
            
            if nonzero_R.size <=100 or nonzero_M.size <= 100 or nonzero_P.size <= 100:
                plt.close(fig_R)
                plt.close(fig_M)
                plt.close(fig_P)
                startHiC = EndHiC
                startSig += Step
                continue
            
            vmaxR = np.percentile(nonzero_R,95)
            vmaxM = np.percentile(nonzero_M,95)
            vmaxP = np.percentile(nonzero_P,95)
            vmaxA = max(vmaxM,vmaxP)
            
            scR = axR1.imshow(Matrix_R, cmap = my_cmap, aspect = 'auto',interpolation = 'none',
                           extent = (0,interval,0,interval),vmax = vmaxR,origin = 'lower')
            scM = axM1.imshow(Matrix_M, cmap = my_cmap, aspect = 'auto',interpolation = 'none',
                           extent = (0,interval,0,interval),vmax = vmaxA,origin = 'lower')
            scP = axP1.imshow(Matrix_P, cmap = my_cmap, aspect = 'auto',interpolation = 'none',
                           extent = (0,interval,0,interval),vmax = vmaxA,origin = 'lower')
            
            for lp in tmp_loops:
                if lp['M_C'] > lp['P_C']:
                    mark = 'M'
                else:
                    mark = 'P'
                x = lp['start'] // res - startHiC
                y = lp['end'] // res - startHiC
                axR1.scatter(y+0.5,x+0.5, color = '', edgecolors = 'b', s = 10)
                axR1.text(y-0.5,x-2.5,'%s' % mark,size = 12)
                axM1.scatter(y+0.5,x+0.5, color = '', edgecolors = 'b', s = 10)
                axM1.text(y-0.5,x-2.5,'%s' % mark,size = 12)
                axP1.scatter(y+0.5,x+0.5, color = '', edgecolors = 'b', s = 10)
                axP1.text(y-0.5,x-2.5,'%s' % mark,size = 12)
            
            ticks = list(np.linspace(0,interval,5).astype(int))
            pos = [((startHiC + t) * ResHiC) for t in ticks]
            labels = [properU(p) for p in pos]
            for ax in [axR1,axM1,axP1]:
                ax.set_xticks(ticks)
                ax.set_xticklabels(labels)
                ax.set_yticks(ticks)
                ax.set_yticklabels(labels)  
            ax = fig_R.add_axes([Left + width + 0.02, HB, 0.01, HH])
            fig_R.colorbar(scR,cax = ax)
            ax = fig_M.add_axes([Left + width + 0.02, HB, 0.01, HH])
            fig_M.colorbar(scM,cax = ax)
            ax = fig_P.add_axes([Left + width + 0.02, HB, 0.01, HH])
            fig_P.colorbar(scP,cax = ax)
            
            
            #CTCF Sigs
            color = '#E41A1C'
            
            axR2 = fig_R.add_axes([Left,SB,width,SH])
            sigs = CTCF_Sigs[i][startSig:EndSig]
            axR2.fill_between(np.arange(len(sigs)),sigs,color = color)
            Sax_Set(axR2,color)
            axR2.set_xlim(0,(EndSig - startSig))
            axR2.set_ylim(0,sigs.max() * 0.8)
            axR2.set_ylabel('CTCF',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)
            axM2 = fig_M.add_axes([Left,SB,width,SH])
            M_sigs = CTCF_Sigs_M[i][startSig:EndSig]
            P_sigs = CTCF_Sigs_P[i][startSig:EndSig]
            s_max = max(M_sigs.max(),P_sigs.max()) * 0.8
            axM2.fill_between(np.arange(len(M_sigs)),M_sigs,color = color)
            Sax_Set(axM2,color)
            axM2.set_xlim(0,(EndSig - startSig))
            axM2.set_ylim(0,s_max)
            axM2.set_ylabel('CTCF',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)            
            axP2 = fig_P.add_axes([Left,SB,width,SH])
            axP2.fill_between(np.arange(len(P_sigs)),P_sigs,color = color)
            Sax_Set(axP2,color)
            axP2.set_xlim(0,(EndSig - startSig))
            axP2.set_ylim(0,s_max)
            axP2.set_ylabel('CTCF',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)
            
            SB += SH
            
            #ATAC Sigs
            color = '#377EB8'

            axR2 = fig_R.add_axes([Left,SB,width,SH])
            sigs = ATAC_Sigs[i][startSig:EndSig]
            axR2.fill_between(np.arange(len(sigs)),sigs,color = color)
            Sax_Set(axR2,color)
            axR2.set_xlim(0,(EndSig - startSig))
            axR2.set_ylim(0,sigs.max() * 0.8)
            axR2.set_ylabel('DNase',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)
            axM2 = fig_M.add_axes([Left,SB,width,SH])
            M_sigs = ATAC_Sigs_M[i][startSig:EndSig]
            P_sigs = ATAC_Sigs_P[i][startSig:EndSig]
            s_max = max(M_sigs.max(),P_sigs.max()) * 0.8
            axM2.fill_between(np.arange(len(M_sigs)),M_sigs,color = color)
            Sax_Set(axM2,color)
            axM2.set_xlim(0,(EndSig - startSig))
            axM2.set_ylim(0,s_max)
            axM2.set_ylabel('DNase',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)            
            axP2 = fig_P.add_axes([Left,SB,width,SH])
            axP2.fill_between(np.arange(len(P_sigs)),P_sigs,color = color)
            Sax_Set(axP2,color)
            axP2.set_xlim(0,(EndSig - startSig))
            axP2.set_ylim(0,s_max)
            axP2.set_ylabel('DNase',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)          
            
            SB += SH
            
            #RNA Sigs
            color = '#4DAF4A'
            
            axR2 = fig_R.add_axes([Left,SB,width,SH])
            sigs = ATAC_Sigs[i][startSig:EndSig]
            axR2.fill_between(np.arange(len(sigs)),sigs,color = color)
            Sax_Set(axR2,color)
            axR2.set_xlim(0,(EndSig - startSig))
            axR2.set_ylim(0,sigs.max() * 0.8)
            axR2.set_ylabel('RNA',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)
            axM2 = fig_M.add_axes([Left,SB,width,SH])
            M_sigs = RNA_Sigs_M[i][startSig:EndSig]
            P_sigs = RNA_Sigs_P[i][startSig:EndSig]
            s_max = max(M_sigs.max(),P_sigs.max()) * 0.8
            axM2.fill_between(np.arange(len(M_sigs)),M_sigs,color = color)
            Sax_Set(axM2,color)
            axM2.set_xlim(0,(EndSig - startSig))
            axM2.set_ylim(0,s_max)
            axM2.set_ylabel('RNA',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)            
            axP2 = fig_P.add_axes([Left,SB,width,SH])
            axP2.fill_between(np.arange(len(P_sigs)),P_sigs,color = color)
            Sax_Set(axP2,color)
            axP2.set_xlim(0,(EndSig - startSig))
            axP2.set_ylim(0,s_max)
            axP2.set_ylabel('RNA',labelpad = 50,rotation = 'horizontal',
                            style = 'italic',size = 10)           
            
            startHiC = EndHiC
            startSig = EndSig
            
            SB = HB + HH
            pp_R.savefig(fig_R)
            pp_M.savefig(fig_M)
            pp_P.savefig(fig_P)
            
            plt.close(fig_R)
            plt.close(fig_M)
            plt.close(fig_P)
            
        pp_R.close()
        pp_M.close()
        pp_P.close()
        

    
