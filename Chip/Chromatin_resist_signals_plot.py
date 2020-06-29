# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 19:02:40 2020

@author: han-luo
"""

def pca_To_200K(fil):
    """
    """
    pca_type = np.dtype({'names':['chr' , 'PCA'],
                     'formats':['S4' , np.float]})
    PCA_Data = np.loadtxt(fil , dtype=pca_type)
    
    chroms = set(PCA_Data['chr'])
    New_Data = {}
    for c in chroms:
        New_Data[c] = {}
        tmp_data = PCA_Data[PCA_Data['chr'] == c]
        New_Data[c] = []
        for i in tmp_data:
            New_Data[c].extend([i['PCA']])
            
    
    return New_Data

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


def Get_PCA_ax(fig , n , PCA , color , labels , ylim):
    ##PCA Tracks
    PCA_index , PCA = UpdateDI(PCA)
    ax = fig.add_axes([Left  , HB + 0.07 * n , width , 0.07])
    ax.fill_between(PCA_index , PCA , where = PCA >= 0 , facecolor = '#E47833' , edgecolor = 'none' )
    ax.fill_between(PCA_index , PCA , where = PCA <= 0 , facecolor = '#7093DB' , edgecolor = 'none' )
    ax.set_xlim(0 , PCA_index.max())
    ax.set_ylim(-0.05 , ylim)
    ax.set_ylabel(labels + '_PC1')
    caxis_S_horizontal(ax,color)

def Get_signal_ax(fig , n , sig , color , labels , ylim):
    ax = fig.add_axes([Left  , HB + 0.07 * n , width , 0.07])
    ax.fill_between(np.arange(len(sig)),sig, facecolor = color, edgecolor = 'none')
    ax.set_xlim((0 , len(sig)))
    ax.set_ylim((0 , ylim))
    ax.set_ylabel(labels,fontsize = 15,rotation = 'horizontal',labelpad = 50)
    caxis_S_horizontal(ax,color)
    return ax
    
size = (12, 12)
Left = 0.2 ; HB = 0.15 ; width = 0.7 ; HH = 0.7

    
CCS_PCA = pca_To_200K('/public/home/xxli/data/BDF1_New/HiC/Compartment/CCS_compartment_200K.txt')
NT5_PCA = pca_To_200K('/public/home/xxli/data/BDF1_New/HiC/Compartment/NT5_compartment_200K.txt')
NT6_PCA = pca_To_200K('/public/home/xxli/data/BDF1_New/HiC/Compartment/NT6_compartment_200K.txt')
fESC_PCA = pca_To_200K('/public/home/xxli/data/BDF1_New/HiC/Compartment/fESC_compartment_200K.txt')

chro = '7'
ccs_pca = np.array(CCS_PCA[chro][23600000 // 200000:31600000 // 200000])
nt5_pca = np.array(NT5_PCA[chro][23600000 // 200000:31600000 // 200000])
nt6_pca = np.array(NT6_PCA[chro][23600000 // 200000:31600000 // 200000])
fesc_pca = np.array(fESC_PCA[chro][23600000 // 200000:31600000 // 200000])


pp = PdfPages('/public/home/xxli/data/BDF1_New/HiC/Compartment/compartment_classify_byself_new/H3K9me3_marked_compartment/Compartment_c7_H3K9me3_RNA_ATAC_signals_26.6_26.8.pdf')
start = 26600000 // 100
end = 26800000 // 100
sig1 = chips[chro][start:end]
sig2 = inputs[chro][start:end]
sig3 = ccs_rna[chro][start:end]
sig4 = nt5_rna[chro][start:end]
sig5 = nt6_rna[chro][start:end]
sig6 = fesc_rna[chro][start:end]
sig7 = ccs_pca
sig8 = nt5_pca
sig9 = nt6_pca
sig10 = fesc_pca
rna_max = max([sig3.max() , sig4.max() , sig5.max() , sig6.max()]) + 5
atac_max = max([sig7.max() , sig8.max() , sig9.max() , sig10.max()]) + 5
fig = plt.figure(figsize = size)
ax = fig.add_axes([Left  , HB , width , 0])
xtick = [0 , len(sig1)]
labels = [properU(start * 100) , properU(end * 100)]
ax.set_xticks(xtick)
ax.set_xticklabels(labels)
ax.set_xlabel('chr' + chro)
ax1 = Get_signal_ax(fig , 0 , sig2 , 'green' , 'K9_input' , sig1.max())
ax2 = Get_signal_ax(fig , 1 , sig1 , 'green' , 'K9_chip' , sig1.max())
ax3 = Get_signal_ax(fig , 2 , sig6 , 'crimson' , 'fESC_RNA' , rna_max)
ax4 = Get_signal_ax(fig , 3 , sig5 , 'crimson' , 'NT6_RNA' , rna_max)
ax5 = Get_signal_ax(fig , 4 , sig4 , 'crimson' , 'NT5_RNA' , rna_max)
ax6 = Get_signal_ax(fig , 5 , sig3 , 'crimson' , 'CCS_RNA' , rna_max)
ax7 = Get_PCA_ax(fig , 6 , sig10 , 'blue' , 'fESC_ATAC' , 0.08)
ax8 = Get_PCA_ax(fig , 7 , sig9 , 'blue' , 'NT6_ATAC' , 0.08)
ax9 = Get_PCA_ax(fig , 8 , sig8 , 'blue' , 'NT5_ATAC' , 0.08)
ax10 = Get_PCA_ax(fig , 9 , sig7 , 'blue' , 'CCS_ATAC' , 0.08)
pp.savefig(fig)
pp.close()
    
    

    
