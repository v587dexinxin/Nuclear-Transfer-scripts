data1 <- cbind(floor(apply(data[c('CCS_R1' , 'CCS_R2' , 'CCS_R3')], 1, mean)) , floor(apply(data[c('NT5_R1' , 'NT5_R2' , 'NT5_R3')], 1, mean)) , 
               floor(apply(data[c('NT6_R1' , 'NT6_R2' , 'NT6_R3')], 1, mean)) , floor(apply(data[c('F35_R1' , 'F35_R2' , 'F35_R3')], 1, mean)) , 
               floor(apply(data[c('F40_R1' , 'F40_R2' , 'F40_R3')], 1, mean)))

colnames(data1) <- c('CCS' , 'NT5' , 'NT6' , 'F35' , 'F40')
sample1 <- sample[c('CCS_R1' , 'NT5_R1' , 'NT6_R1' , 'F35_R1' , 'F40_R1'),]
row.names(sample1) <- c('CCS' , 'NT5' , 'NT6' , 'F35' , 'F40')
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = data1, colData = sample1,  design= ~ batch + conditions)
dds <- DESeq(ddsFullCountTable)


normalized=T
normalized_counts <- counts(dds, normalized=TRUE)
head(normalized_counts)
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
normalized_counts



rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]


hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(rlogMat, method="pearson"))
hc <- hcluster(t(rlogMat), method="pearson")
hM <- format(round(pearson_cor, 3))


heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none", col=hmcol, margins=c(11,11), main="The pearson correlation of eachsample" , cellnote = hM,notecol='black')


##---------------------remove_batch effects----------------------



vsd <- vst(dds, blind=FALSE)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)


hmcol <- colorRampPalette(c("blue","white","red"))
pearson_cor <- as.matrix(cor(assay(vsd), method="pearson"))
hc <- hcluster(t(assay(vsd)), method="pearson")
hM <- format(round(pearson_cor, 3))


pdf("F:/work/ntESC_3Dreprogramming/Workspace_New/data/RNA/RNA_New/Correlation_plot/RNA_correlation_heatmap_remove_batch.pdf", pointsize=10)
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none", col=hmcol, margins=c(11,11), main="The pearson correlation of eachsample" , cellnote = hM,notecol='black')
dev.off()
















ggplot(pcaData, aes(PC1, PC2, color=conditions, shape = batch))
geom_point(size=3)
xlim(-12, 12)
ylim(-10, 10)
xlab(paste0("PC1: ",percentVar[1],"% variance"))
ylab(paste0("PC2: ",percentVar[2],"% variance"))
geom_text(aes(label=name),vjust=2)



















