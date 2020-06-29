library(splitstackshape)
gene.id <- read.table(file.choose())
gene.id <- cSplit(gene.id,"V1",sep = ".", drop = F)

gene.id.char <- as.character(gene.id$V1_1)

A <- gene.id$V1_1

fESC_data = read.table('D:\\Workspace_New\\data\\ATAC\\peaks\\fESC_R1_R2_Common_peaks.txt')
fESC_data = fESC_data[,8:9]
fESC_M = as.matrix(fESC_data)
fESC_out = est.IDR(fESC_M, mu, sigma, rho, p, eps=0.001,max.ite = 20)
fESC.selected = select.IDR(fESC_M,fESC_out$IDR,0.01)
write.table(fESC.selected$x,'D:\\Workspace_New\\data\\ATAC\\peaks\\fESC_R1_R2_selected.txt',row.names = FALSE,quote = FALSE,col.names = FALSE)