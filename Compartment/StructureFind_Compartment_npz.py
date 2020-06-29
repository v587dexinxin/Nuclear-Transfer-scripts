# -*- coding: utf-8 -*-
"""
Created on Fri Nov 01 16:20:19 2019

@author: han-luo
"""

from __future__ import division
import matplotlib
matplotlib.use('Agg')
from scipy import sparse
from sklearn import  isotonic
from sklearn.decomposition import PCA
from collections import OrderedDict
from scipy.stats import poisson
from statsmodels.sandbox.stats.multicomp import multipletests
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from itertools import  islice
import os, math, random, ghmm, bisect
import cooler, copy
import numpy as np


# A Customized  Chromotin Structure Finding Class
class StructureFind(object):
    """
    This is  a module class for Chromosome 3D Structure Analysis based on Hi-C Data.
    The basis class can fix the haplotype and non-haplotype structures.
    We can use this class to Get :
    
                                  Compartments
                                  TADs.
                                  Chromatin Loops.
        
    This class includes many self-defined methods. But if U do not want fussy details,
    Only some key Functions API are needed.
    
    Input Data
    ----------
    
    This module follows the Matrix building by HiCHap sub-command (matrix).
    The results of Interaction Matrix will be saved as cooler format.
        
    
    
    ---------------------------------------------------------------------------------------------
    
    More API information can be found in methods.
    
    code::block 
    
    	>>> from HiCHap.StructureFind import StructureFind
	
	>>> #============= Compartment==============
	>>> ## For traditional Hi-C
	>>> GM_T_PC = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 500000, Allelic = False)
	>>> GM_T_PC.run_Compartment(OutPath = 'Traditonal_PC', plot = True, MS = 'IF', SA = False)
	 
	>>> ## For haplotype-resolved Hi-C 
	>>> GM_M_PC = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 500000, Allelic = 'Maternal')
	>>> GM_M_PC.run_Compartment(OutPath = 'Maternal_PC', plot = True, MS = 'IF', SA = False)
	
	>>> GM_P_PC = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 500000, Allelic = 'Paternal')
	>>> GM_P_PC.run_Compartment(OutPath = 'Paternal_PC', plot = True, MS = 'IF', SA = False)

	
	>>> #============= TADs calling=============
	>>> ## For traditional Hi-C 
	>>> GM_tads_T = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 40000, Allelic = False)
	>>> GM_tads_T.run_TADs(OutPath = 'Traditional_TADs', plot = True)
	>>>
	>>> ## For haplotype-resolved Hi-C
	>>> GM_tads_M = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Maternal')
	>>> GM_tads_M.run_TADs(OutPath = 'Maternal_TADs', plot = True)

	>>> GM_tads_P = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Paternal')
	>>> GM_tads_P.run_TADs(OutPath = 'Paternal_TADs', plot = True)
	

	>>> #============= Loops calling=============
	>>> ## For traditonal Hi-C
	>>> GM_Loop_T = StructureFind(cooler_fil = 'Merged_Traditional_Multi.cool', Res = 40000, Allelic = False)
	>>> GM_Loop_T.run_Loops(OutPath = 'Traditional_Loops', plot = True)
	
	>>> ## For haplotype-resolved Hi-C
	>>> GM_Loop_M = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Maternal')
	>>> GM_Loop_M.run_Loops(OutPath = 'Maternal_Loops', plot = True)
	
	>>> GM_Loop_P = StructureFind(cooler_fil = 'Merged_Imputated_Haplotype_Multi.cool', Res = 40000, Allelic = 'Paternal')
	>>> GM_Loop_P.run_Loops(OutPath = 'Paternal_Loops', plot = True)
	
 
    
    """
    def __init__(self, npz_fil, Res, Allelic, GapFile = None, 
                 Loop_ratio = 0.6, Loop_strength = 16):
        
        #------------Initialization the HiC Matrix Dictionary----------------
        self.npz_fil = npz_fil
        self.Res = Res
        self.Allelic = Allelic
        self.Gap_file = GapFile
        self.ratio = Loop_ratio
        self.LoopStrength = Loop_strength

#===============================Public Functions====================================


    def Sig_update(self, sigs):
        """
            Up date the Sigs to Plot.
            
        Parameter
        ---------
        sigs : array
            sigs array.
        """
        New_Sig = []
        New_Index = []
        
        for index in range(len(sigs) - 1):
            if sigs[index] * sigs[index + 1] < 0:
                New_Sig.append(sigs[index])
                New_Sig.append(0)
                New_Index.append(index)
                New_Index.append(index + 0.5)
            else:
                New_Sig.append(sigs[index])
                New_Index.append(index)
        
        return np.array(New_Index), np.array(New_Sig)
    
    
    def Cmap_Setting(self, types = 2, **kwargs):
        """
            Color bar setting. hexadecimal notation(#16) must be input.
        
        Parameters
        ----------
        start_Color : str
            start_Color (default : #FFFFFF)        
        
        end_Color : str
            end_Color (default : #CD0000)
        """
        
        start_Color = kwargs.get('start_Color','darkblue')
        middle_Color = kwargs.get('middle_Color','black')
        end_Color = kwargs.get('end_Color','red')
        if types == 2:
            return LinearSegmentedColormap.from_list('interactions',[start_Color,end_Color])
        elif types == 3:
            return LinearSegmentedColormap.from_list('interactions',[start_Color,middle_Color,end_Color])
            


    def properU(self, pos):
        """
          Express a genomic position in a proper unit (KB, MB, or both).
          
        """
        i_part = int(pos) // 1000000    # Integer Part
        d_part = (int(pos) % 1000000) // 1000   # Decimal Part
        
        if (i_part > 0) and (d_part > 0):
            return ''.join([str(i_part), 'M', str(d_part), 'K'])
        elif (i_part == 0):
            return ''.join([str(d_part), 'K'])
        else:
            return ''.join([str(i_part), 'M'])
            
            
    def caxi_H(self,ax):
        """
            Axis Control for HeatMaps
            
        """
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        ax.tick_params(axis = 'both', labelsize = 10, length = 5, pad = 7)
    
    def fig_settings(self):
        """
            About figure.
        """
        self.figure_Dict = {}
        self.figure_Dict['size'] = (10,9)
        self.figure_Dict['width'] = 0.618
        self.figure_Dict['Left'] = (1 - self.figure_Dict['width']) / 2
        self.figure_Dict['HB'] = 0.1
        self.figure_Dict['HH'] = self.figure_Dict['width'] * self.figure_Dict['size'][0] / self.figure_Dict['size'][1]
        self.figure_Dict['SB'] = self.figure_Dict['HB'] + self.figure_Dict['HH']
        
    
#============================Compartment Analysis modules===========================


    
    def Distance_Decay(self, M, G_array):
        """
            Get the distance decay line.
            
        Parameter
        ---------
        M : 2D HiC Matrix of Intra-Chromosome
        
        G_array : Gap if not available, Get it.
        
        """
        
        size = M.shape[0]
        bin_arange = np.arange(size)        
        
        if G_array == None:
            Gap_ratio = 0.05
            nonzero_mask = (M != 0).sum(axis = 0) / float(size)
            gap_mask = np.where(nonzero_mask <= Gap_ratio,True,False)
            G_array = bin_arange[gap_mask]
            NG_array = bin_arange[~gap_mask]
        else:
            NG_array = np.array([item for item in bin_arange if item not in G_array])
            
        Matrix_mask = []
        dtype = [('bin1',np.int),('bin2',np.int),('IF',np.float)]
        bin1, bin2 = np.nonzero(M)
        IF = M[bin1,bin2]
        
        Sparse_M = np.zeros((len(IF),),dtype = dtype)
        Sparse_M['bin1'] = bin1
        Sparse_M['bin2'] = bin2
        Sparse_M['IF'] = IF
        
        for item in G_array:
            Matrix_mask.extend(np.where(Sparse_M['bin2'] == item)[0])
        
        Matrix_mask_set = np.array(list(set(Matrix_mask)))
        Mask_True = np.ones((len(Sparse_M),),dtype = bool)
        
        if Matrix_mask_set.shape[0] >= 1:
            Mask_True[Matrix_mask_set] = False
        
        Sparse_M_Less = Sparse_M[Mask_True]
        
        weight = np.array(Sparse_M_Less['IF'])
        weight = np.hstack((weight,np.array([0])))

        distance = np.array(Sparse_M_Less['bin2'] - Sparse_M_Less['bin1'], dtype = np.int)
        distance_abs = np.array([abs(x) for x in distance], dtype = np.int)
        distance_size_abs = np.hstack((distance_abs,np.array([size],dtype = np.int)))
        distance_bin = np.bincount(distance_size_abs,weight)
        
        #-------------------------Remove the Gap effect-----------------------------
        for i in range(size):
            if i == 0:
                gap_num_start = 0
                gap_num_end = sum((0 <= G_array) & (G_array <= size - 1))
                gap_num = gap_num_start + gap_num_end
                bin_num = float(size - i) - gap_num
            else:
                gap_num_start = sum((0 <= G_array) & (G_array <= size - 1 - i))
                gap_num_end = sum((i <= G_array) & (G_array <= size - 1))
                gap_num = gap_num_start + gap_num_end
                bin_num = float(size - i) * 2 - gap_num
            
            if bin_num > 0:
                distance_bin[i] = float(distance_bin[i] / bin_num)
        distance_bin = distance_bin[: size]
    
        return distance_bin, G_array, NG_array
        
    
    def Sliding_Approach(self, M, decline, window):
        """
            Sliding Approach by Ren Bin.  may be a little difference

        """
        step = window // self.Res // 2
        
        N = M.shape[0]
        
        OE_matrix_bigger = np.zeros(M.shape)
        
        for i in range(N):
            for j in range(N):
                if i < step or j < step:
                    OE_matrix_bigger[i][j] = M[i][j] / decline[abs(i-j)]
                elif i > N - step - 1 or j > N - step - 1:
                    OE_matrix_bigger[i][j] = M[i][j] / decline[abs(i-j)]
                else:
                    O_Sum = M[i-step:i+step+1, j-step:j+step+1].sum()
                    E_Sum = 3 * decline[abs(i-j)] + 2 * decline[abs(i-j-1)] + \
                            2 *decline[abs(i-j+1)] + decline[abs(i-j-2)] + \
                            decline[abs(i-j+2)]
                    
                    OE_matrix_bigger[i][j] = O_Sum / E_Sum
        
        return OE_matrix_bigger
                    
        
    def Get_PCA(self, distance_bin, M, NG_array, SA = False):
        """
            Based on Distance_Decay 
            Get the:
                Correlation Matrix
                Obersevered / Expected Matrix (OE Matrix)
                PCA component
        
        Parameters
        ----------
        distance_bin : distance_bin created by Distance_Decay
        
            M        : 2D HiC Matrix of Intra-Chromosome
            
        NG_array     : NG_array created by Distance_Decay
        
        """
        
        decline = distance_bin
        decline[decline == 0] = decline[np.nonzero(decline)].min()
        
        OE_matrix_bigger = np.zeros(M.shape)
        
        if SA == False:
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    if M[i][j] != 0:
                        OE_matrix_bigger[i][j] = M[i][j] / decline[abs(i-j)]
        else:
            OE_matrix_bigger = self.Sliding_Approach(M, decline, window=600000)
        
        OE_matrix = OE_matrix_bigger[:,NG_array]
        
        Cor_matrix = np.corrcoef(OE_matrix,rowvar = False)
        Cor_matrix[np.isnan(Cor_matrix)] = 0
        Cor_matrix[np.isinf(Cor_matrix)] = 1
        pca = PCA(n_components = 3)
        pca.fit(Cor_matrix)
        pca_components=np.vstack((pca.components_[0],pca.components_[1],pca.components_[2]))
        
        return pca_components, Cor_matrix, OE_matrix
        
    
    def Select_PC(self, Cor_matrix, pca_components):
        """
            Based On Get_PCA
            Select the PC 
        
        Parameters
        ----------
        Cor_matrix : Cor_matrix created by Get_PCA
        
        OE_matrix  : OE_matrix created by Get_PCA
        
        """
        select_k = 0
        corr_coef = 0
        for i in range(pca_components.shape[0]):
            tmp_coef = np.array([np.corrcoef(pca_components[i], item)[0,1] for item in Cor_matrix])
            tmp_coef[np.isnan(tmp_coef)] = 0
            tmp_coef[np.isinf(tmp_coef)] = 1
            if abs(tmp_coef).sum() > corr_coef:
                corr_coef = abs(tmp_coef).sum()
                select_k = i
                direction = 1
                if tmp_coef.sum() < 0:
                    direction = -1
        
        pc_selected = pca_components[select_k] * direction
        
        return pc_selected
    
    
    def Loading_Tranditional_PC(self, fil):
        """
        """
        PC_Dic = {}
        with open(fil,'r') as f:
            for line in f:
                line = line.strip().split()
                if line[0] not in PC_Dic.keys():
                    PC_Dic[line[0]] = []
                    PC_Dic[line[0]].append(line[-1])
                else:
                    PC_Dic[line[0]].append(line[-1])
        
        for chro, v in PC_Dic.items():
            v = np.array(v, dtype = float)
            PC_Dic[chro] = v
        
        return PC_Dic
    
    
    def Select_Allelic_PC(self, pca_components, Tranditional_PC, eps = 0.7):
        """
            Select the principal component by  supervised manner in Allelic matrix
            
        """
        PCC = []
        for pc in pca_components:
            pcc = abs(np.corrcoef(pc, Tranditional_PC)[0][1])
            PCC.append(pcc)
        if np.max(PCC) < eps:
            print "    PCC too low for this chromosome, check it if possible!"

        index = np.argmax(PCC)
    
        return pca_components[index]
     
     
    def Refill_Gap(self, M1, M2, NonGap, dtype):
        """
            Refill zero value into Gap in OE Matrix and Correlation Matrix.
            
        """
        Refilled_M = np.zeros(M1.shape, dtype = np.float)
        if dtype == 'Cor':
            tmp = np.zeros((M1.shape[0], M2.shape[0]), dtype = np.float)
            for i in range(len(NonGap)):
                insert_pos = NonGap[i]
                tmp[insert_pos] = M2[i]

            tmp = tmp.T
        
            for i in range(len(NonGap)):
                insert_pos = NonGap[i]
                Refilled_M[insert_pos] = tmp[i]

        elif dtype == 'OE':   
            M2 = M2.T
            for i in range(len(NonGap)):
                insert_pos = NonGap[i]
                Refilled_M[insert_pos] = M2[i]
                
                Refilled_M = Refilled_M.T
                
        return Refilled_M
        
    def Compartment(self, SA = False, Tranditional_PC_file = None):
        """
            Compute the Compartment Structure for each Chromosome
        """
        
        print "Compartment Analysis start ..."
        print "Resolution is %s ..." % self.properU(self.Res)
        
        self.npz = np.load(self.npz_fil , allow_pickle=True)
        self.npz = self.npz['Matrix'][()]
        
        if self.Allelic == False:
            chroms = self.npz.keys()
        elif self.Allelic == 'Maternal':
            chroms = [i for i in self.npz.keys() if i.startswith('M')]
        elif self.Allelic == 'Paternal':
            chroms = [i for i in self.npz.keys() if i.startswith('P')]
        else:
            raise Exception ('Unkonwn key word %s, Only Maternal, Paternal, False allowed' % self.Allelic)
        
        self.chroms = chroms
        Matrix_Lib = {}
        for chro in chroms:
            Matrix_Lib[chro] = self.npz[chro]
        
        self.Matrix_Dict = Matrix_Lib
        self.Cor_Martrix_Dict = {}
        self.OE_Matrix_Dict = {}
        self.Compartment_Dict = {}
        
        if self.Allelic != False:
            Tranditional_PC = self.Loading_Tranditional_PC(Tranditional_PC_file)
        
        self.RawPCA = {}
        for chro in chroms:
            print "Chromosome %s start" % chro
            M = Matrix_Lib[chro]
            self.RawPCA[chro] = []
            distance_bin, Gap, NonGap = self.Distance_Decay(M=M, G_array=None)
            pca, Cor_M, OE_M = self.Get_PCA(distance_bin=distance_bin, M=M, NG_array=NonGap, SA = SA)
            
            if self.Allelic == False:
                pc_select = self.Select_PC(Cor_matrix=Cor_M, pca_components=pca)
                Compartment_zeros = np.zeros((M.shape[0],), dtype = np.float)
                Compartment_zeros[NonGap] = pc_select
                self.Compartment_Dict[chro] = Compartment_zeros
                
            else:
                for i in range(len(pca)):
                    tmp = np.zeros((M.shape[0],), dtype = np.float)
                    tmp[NonGap] = pca[i]
                    self.RawPCA[chro].append(tmp.tolist())

                self.RawPCA[chro] = np.array(self.RawPCA[chro])
                Compartment_zeros = np.zeros((M.shape[0],), dtype = np.float)
                tran_chro = chro[1:]
                pc_select = self.Select_Allelic_PC(self.RawPCA[chro], Tranditional_PC[tran_chro])
                Compartment_zeros[NonGap] = pc_select[NonGap]
                self.Compartment_Dict[chro] = Compartment_zeros
                
            OE_M = self.Refill_Gap(M, OE_M, NonGap, dtype = 'OE')
            Cor_M = self.Refill_Gap(M, Cor_M, NonGap, dtype = 'Cor')

            self.Cor_Martrix_Dict[chro] = Cor_M
            self.OE_Matrix_Dict[chro] = OE_M
   

    def OutPut_PC_To_txt(self, out):
        """
            Output the PC component to a text file
        
        Parameter
        ---------
        out : str
            out file
        """
        f = open(out, 'w')
        if self.Allelic == False:
            for chro, pc in self.Compartment_Dict.items():
                for value in pc:
                    f.writelines(chro+'\t'+str(value)+'\n')
        else:
            for chro, pc in self.Compartment_Dict.items():
                for value in pc:
                    f.writelines(chro[1:]+'\t'+str(value)+'\n')
        
        f.close()
    
    
    def Plot_Compartment(self, out, MS = 'IF'):
        """
            Output the Compartment Structure to a  PDF
        
        Parameter
        --------
        out : str
            out PDF file.
        
        Matrix :  Matrix Str type
        
                'IF' : Interaction Matrix.
                'OE' : Observed / Expected Matrix.
                'Cor' : Correlation Matrix.
                
        """
        pp = PdfPages(out)
        self.fig_settings()
        if MS == 'IF':
            Matrix_Dict = self.Matrix_Dict
            cmap = self.Cmap_Setting()
        elif MS == 'OE':
            Matrix_Dict = self.OE_Matrix_Dict
            cmap = self.Cmap_Setting(types=3,start_Color='darkblue')
        elif MS == 'Cor':
            Matrix_Dict = self.Cor_Martrix_Dict
            cmap = self.Cmap_Setting(types=3,start_Color='darkblue')
        else:
            raise Exception("Unknown Matrix String %s, Only IF, OE, Cor strings allowed." % MS)
        
        for chro in self.chroms:
            Matrix = Matrix_Dict[chro]
            Sigs = self.Compartment_Dict[chro]
            N = Matrix.shape[0]
            
            if MS == 'IF':
                nonzero = Matrix[np.nonzero(Matrix)]
                vmax = np.percentile(nonzero,95)
                vmin = 0
            elif MS == 'OE':
                nonzero = Matrix[np.nonzero(Matrix)]
                vmax = np.percentile(nonzero, 90)
                vmin = 2 - vmax
            elif MS == 'Cor':
                nonzero = Matrix[np.nonzero(Matrix)]
                vmax = np.percentile(nonzero, 90)
                vmin = -vmax
            
            fig = plt.figure(figsize = self.figure_Dict['size'])
            ax = fig.add_axes([self.figure_Dict['Left'], self.figure_Dict['HB'],
                               self.figure_Dict['width'], self.figure_Dict['HH']])
                               
            sc = ax.imshow(Matrix, cmap = cmap, aspect = 'auto', interpolation = 'none',
                           extent = (0, N, 0, N), vmin = vmin, vmax = vmax, origin = 'lower')
            
            ticks = list(np.linspace(0, N, 5).astype(int))
            pos = [t * self.Res for t in ticks]
            labels = [self.properU(p) for p in pos]
            
            ax.set_xticks(ticks)
            ax.set_xticklabels(labels)
            ax.set_yticks(ticks)
            ax.set_yticklabels(labels)
            if self.Allelic == False:
                ax.set_xlabel('Chr'+chro, size = 14)
            else:
                ax.set_xlabel('Chr'+chro[1:], size = 14)
            
            
            
            ax = fig.add_axes([self.figure_Dict['Left']+self.figure_Dict['width']+0.02,
                               self.figure_Dict['HB'],0.01,self.figure_Dict['HH']])
            
            fig.colorbar(sc, cax = ax)
            
            index, sigs = self.Sig_update(Sigs)
            
            ax2 = fig.add_axes([self.figure_Dict['Left'],self.figure_Dict['SB'],
                                self.figure_Dict['width'],self.figure_Dict['HB']])

            for spine in ['right', 'top', 'left']:
                ax2.spines[spine].set_visible(False)
                
            ax2.fill_between(index, sigs, where = sigs <= 0, color = 'midnightblue')
            ax2.fill_between(index, sigs, where = sigs >= 0, color = 'gold')
            ax2.tick_params(axis = 'both',bottom = False,top = False,left = False,
                    right = False, labelbottom =False, labeltop = False,
                    labelleft = False, labelright = False)
            
            ax2.set_xlim(0,len(Sigs))
            ax2.set_ylabel('PC', size = 12)
            
            pp.savefig(fig)
            plt.close(fig)
        
        pp.close()                    
        
        
    def run_Compartment(self, OutPath, plot = True, MS = 'IF', 
                        SA = False, Tranditional_PC_file = None):
        """
            Main function to Get Compartment
        
        Parameters
        ----------
        OutPath : str
            Out Put Path
        """
        OutPath.rstrip('/')
        if os.path.exists(OutPath):
            pass
        else:
            os.mkdir(OutPath)
        
        prefix = os.path.split(OutPath)[-1]
        res = self.properU(self.Res)
        
        fil = os.path.join(OutPath, prefix+'_Compartment_'+res+'.txt')
        pdf = os.path.join(OutPath, prefix+'_Compartment_'+MS+'_'+res+'.pdf')
        
        self.Compartment(SA = SA, Tranditional_PC_file=Tranditional_PC_file)
        self.OutPut_PC_To_txt(fil)
        if plot:
            self.Plot_Compartment(pdf, MS)