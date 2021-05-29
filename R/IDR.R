library(splitstackshape)
gene.id <- read.table(file.choose())
gene.id <- cSplit(gene.id,"V1",sep = ".", drop = F)

gene.id.char <- as.character(gene.id$V1_1)

A <- gene.id$V1_1


mu <- 2.6
sigma <- 1.3
rho <- 0.8
p <- 0.7

library('idr')
fESC_data = read.table('F:/work/ntESC_3Dreprogramming/Workspace_New/data/ATAC/ATAC_new/peaks/idr_new/NT5_R1_R2_Common_peaks.narrowPeak')
fESC_data <- fESC_data[,c('V8','V18')]
fESC_M = as.matrix(fESC_data)
fESC_out = est.IDR(fESC_M, mu, sigma, rho, p, eps=0.001,max.ite = 20)
fESC.selected = select.IDR(fESC_M,fESC_out$IDR,0.01)
write.table(fESC.selected$x,'F:/work/ntESC_3Dreprogramming/Workspace_New/data/ATAC/ATAC_new/peaks/idr_new/CCS_R1_R2_selected.txt',row.names = FALSE,quote = FALSE,col.names = FALSE)



CCS.selected = select.IDR(CCS_M,fESC_out$IDR,0.01)
NT5.selected = select.IDR(NT5_M,fESC_out$IDR,0.01)
NT6.selected = select.IDR(NT6_M,fESC_out$IDR,0.01)
F35.selected = select.IDR(F35_M,fESC_out$IDR,0.05)
F40.selected = select.IDR(F40_M,fESC_out$IDR,0.001)






##------------------NTs¡¢fESC---------------------------

mu <- 2.6
sigma <- 1.3
rho <- 0.8
p <- 0.7

library('idr')
fESC_data = read.table('F:/work/ntESC_3Dreprogramming/Workspace_New/data/ATAC/ATAC_new/peaks/idr_new/fESC_common_peaks.narrowPeak')
fESC_data <- fESC_data[,c('V8','V28')]
fESC_M = as.matrix(fESC_data)
fESC_out = est.IDR(fESC_M, mu, sigma, rho, p, eps=0.001,max.ite = 20)
fESC.selected = select.IDR(fESC_M,fESC_out$IDR,0.01)
write.table(fESC.selected$x,'F:/work/ntESC_3Dreprogramming/Workspace_New/data/ATAC/ATAC_new/peaks/idr_new/fESC_selected.txt',row.names = FALSE,quote = FALSE,col.names = FALSE)


