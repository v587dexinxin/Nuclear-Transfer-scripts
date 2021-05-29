# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 21:10:41 2020

@author: xxli
"""


from __future__ import  division
from itertools import  islice
from scipy import optimize
import numpy as np
import bisect, cPickle
import pysam, sys, os


    
    
def Unique_Mapping_Filtering(fil,out):
    """
    """
    unfilter_bam = pysam.AlignmentFile(fil,'rb')
    
    filtered_bam = pysam.AlignmentFile(out,'wb', template = unfilter_bam)
    Total = 0
    Unique = 0
    same_chrom = 0
    unmapped = 0
    for read in unfilter_bam:
        Total += 1
        
        if read.reference_name != read.next_reference_name:
            same_chrom += 1
            continue
            
            
        if read.is_unmapped:
            unmapped += 1
            continue
            
        else:
            chro = read.reference_name.lstrip('chr')
            if not (chro.isdigit() or chro == 'X'):
                continue
            if read.has_tag('XS'):
                first = read.get_tag('AS')
                second = read.get_tag('XS')
                if isinstance(second,int):
                    if first > second:
                        filtered_bam.write(read)
                        Unique += 1
                else:
                    filtered_bam.write(read)
                    Unique += 1
                    
                
            else:
                filtered_bam.write(read)
                Unique += 1
    
    print "Total reads number : %s" % Total    
    print "Unique reads number is %s" % Unique
    print "UnSame_chroms reads number is %s" % same_chrom
    print "Unmapped reads number is %s" % unmapped
        
    unfilter_bam.close()
    filtered_bam.close()
    
    
    
if __name__ == '__main__':
    bam_fil = sys.argv[1]
    Unique_fil = 'Unique_'+os.path.split(bam_fil)[1]
    Unique_Mapping_Filtering(bam_fil,Unique_fil)
