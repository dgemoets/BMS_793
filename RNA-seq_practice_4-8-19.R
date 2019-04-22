#R libraries2=already run standardized experiments
#from workflow on https://f1000research.com/articles/4-1070/v2
load("C:/Users/adt0023/Documents/HumanLncRNA2019_noncoding.RData")
#HumanlncRNA2019_noncoding2 is parsed down to results, no SE


library('DESeq2')
library('dplyr')
library('AnnotationDbi')
library('ggplot2')
library('calibrate')
library("PoiClaClu")
library("pheatmap")
library("RColorBrewer")
library("Rsamtools")
library("ggplot2")
library("genefilter")


#Start with raw reads (FASTQ format)
#Align reads to reference genome (into BAM files)
#Gene-guided transcriptome assembly (GTF)
#Condense into summarized experiment to transfer to DESeq for differential expression analysis

#Preview of dataset
colData(se)

##########################Exploratory analysis#############################################

#Creating DESeqDataSet object
dds <- DESeqDataSet(se, design = ~ cell + dex)
#Filtering out rows with no coints or just one count across all samples
dds <- dds[ rowSums(counts(dds)) > 1, ]
#rlog transformation to stabilize variance across mean
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
#sample distances (Euclidean)
sampleDists <- dist( t( assay(rld) ) )
sampleDists
#sample distances visualized
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
#Poisson distance (measures dissimilarity)
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)
#PCA plot
plotPCA(rld, intgroup = c("dex", "cell"))

#Multidimensional scaling
mdsData <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=dex,shape=cell)) + geom_point(size=3)

##########################Differential Expression Analysis################################
#dds experimental design already defined, 
#can run DESeq on dds object to fit return new dds fitted to differential expression pipeline parameters
dds <- DESeq(dds)
#Building results table
(res <- results(dds))
summary(res)
#metadata on meanings of columns
mcols(res, use.names=TRUE)
#p-adjusted
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
#adjusting log fold change threshold to 1
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

##Plotting
#quick plot of counts
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))
#MA-plot
plotMA(res, ylim=c(-11,11))
#MA-plot with log fold change threshold & labeling individual point with lowest adjusted p value
plotMA(resLFC1, ylim=c(-11,11), main="adjusted log fold change threshold MA plot")
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ],{
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
         text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})
#hist of p-values, excluding genes with small counts
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")

#volcano (from http://genomics.as.wvu.edu/bioinfo_scripts.php?fileToGet=bioscripts/rnaseq.html)
par(pch = 16)
with(res, plot(log2FoldChange, -log10(pvalue), main = "Volcano plot"))
with(subset(res, padj < 0.05), points(log2FoldChange, -log10(pvalue), col = "red"))
with(subset(res, abs(log2FoldChange) > 2), points(log2FoldChange, -log10(pvalue),  col = "orange"))
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 2), points(log2FoldChange,  -log10(pvalue), col = "green"))
# Add legend
legend("bottomleft", legend = c("FDR<0.05", "|LFC|>2", "both"), pch = 16, col = c("red", "orange", "green"))
#txtxy package \/
# Label Extra significant points
with(subset(res, padj < 0.05 & abs(log2FoldChange) > 2), textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = 0.5))
# Label all significant
with(subset(res, padj < 0.05), textxy(log2FoldChange, -log10(pvalue), labs = Gene, cex = 0.5))

#Gene clustering
#heatmap 20 genes with highest variance
topVarGenes <- head(order(rowVars(assay(rld)),decreasing=TRUE),20)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell","dex")])
pheatmap(mat, annotation_col=df)

#workflow has code for exporting to CSV, HTML documents as well as for removing batch effects


###########################Differential expression beyond workflow################################################

###Goal=identify genes uniquely expressed in diabetic &/ mito
#cell.mito<-subset(rld, cell="Mito")
#plotPCA(rld, intgroup= c("dex",cell.mito), ntop=30)
#plotPCA(rld, intgroup=subset(rld$c("dex","cell"), cell="mito"))

#topGene <- rownames(res)[which.min(res$padj)]
#plotCounts(dds, gene=topGene, intgroup=c("dex","cell"[["Mito"]]))


##deseq data set->just mito
#se.mito<-subset(se,cell="Mito")
#se.mito<-filter(se, cell=="Mito")
se.mito<-se[cell="Mito"]
dds.mito<-DESeqDataSet(se.mito,design= ~ cell+ dex) #error model matrix not full rank b/c cell= only one value
dds.mito<-DESeqDataSet(se, design = ~ (cell["Mito"]) + dex) #error b/c dex=different length from cell


#heatmap of mito vs cyto (ignoring db vs non-db)
dds.cell<- DESeqDataSet(se, design= ~cell)
rld.cell<-rlog(dds.cell, blind=FALSE)
poisd.cell<- PoissonDistance(t(counts(dds.cell)))
samplePoisDistMatrix.cell <- as.matrix( poisd.cell$dd )
rownames(samplePoisDistMatrix.cell) <- paste( rld.cell$cell)
colnames(samplePoisDistMatrix.cell) <- NULL
colors<- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix.cell,
         clustering_distance_rows=poisd.cell$dd,
         clustering_distance_cols=poisd.cell$dd,
         col=colors, main = "Mito v Cyto")

#heatmap of db vs non-db (ordered)
dds.dex<- DESeqDataSet(se, design= ~dex)
rld.dex<-rlog(dds.dex, blind=FALSE)
#rld.dex<-rld.dex[order(dex),]
rld.dex.mito<-sub
rld.dex<-rld.dex[order(rld.dex$dex),]
poisd.dex<- PoissonDistance(t(counts(dds.dex)))
samplePoisDistMatrix.dex<- as.matrix( poisd.dex$dd )
rownames(samplePoisDistMatrix.dex) <- paste( rld.dex$dex)
colnames(samplePoisDistMatrix.dex) <- NULL
colors<- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(samplePoisDistMatrix.dex,
         clustering_distance_rows=poisd.dex$dd,
         clustering_distance_cols=poisd.dex$dd,
         col=colors, main = "Db v non-Db")

#heatmap of just mito with db and non-db (?)
sampleDists.mito <- dist( t( assay(subset(rld, cell="Mito"))))
sampleDists.mito
sampleDistMatrix.mito <- as.matrix( sampleDists.mito)
rownames(sampleDistMatrix.mito) <- paste( rld$dex, rld$cell, sep="-" )
colnames(sampleDistMatrix.mito) <- NULL
colors.mito <- colorRampPalette( rev(brewer.pal(9, "Greens")) )(255)
pheatmap(sampleDistMatrix.mito,
         clustering_distance_rows=sampleDists.mito,
         clustering_distance_cols=sampleDists.mito,
         col=colors.mito, main="Mito?")
#mislabelled axis?


##differential expression of a single gene (from book)
#variance stabilized data to simplify considerations
vsd<-getVarianceStabilizedData(dds)
heatmap(cor(vsd),cexCol=0.75,cexRow=0.75, col=colors)

#expression barplot (head)
sig.ordered<- res.05[order(res.05$padj),]
for(gene in head(rownames(sig.ordered))){
  barplot(vsd[gene,],las=2,col=as.numeric(as.factor(gene)),main=gene,cex.names = 0.5)
}
#Expression barplot (tail, least variable gene just for proof of plot)
for(gene in tail(rownames(sig.ordered))){
  barplot(vsd[gene,],las=2,col=as.numeric(as.factor(gene)),main=gene,cex.names = 0.5)
}

#"ENSG00000251562", MALAT1 (gene id=not technically part of a column)
res["ENSG00000251562",]
#mal.df<-as.data.frame(vsd["ENSG00000251562",])  class S4, can't be coerced into df

barplot(vsd["ENSG00000251562",],las=2, main="ENSG00000251562", cex.names=0.5)

for(gene in rownames(vsd)=="ENSG00000251562"){
  barplot(res[gene,],las=2,col=as.numeric(as.factor(gene)),main=gene,cex.names = 0.5)
}



#heatmap specific to MALAT1 (and NEAT1?)
#mat <- assay(rld)[ topVarGenes, ]
#mat <- assay(rld)[rld["ENSG00000251562",]]
mat <- assay(rld)[c("ENSG00000251562", "ENSG00000245532")]
mat <- mat - rowMeans(mat) #x must be an array of at least two dimensions
df <- as.data.frame(colData(rld)[,c("cell","dex")])
pheatmap(mat, annotation_col=df, main="MALAT1 and NEAT1")


##Venn diagram for overlap between treatments and genes
#venn package (not Venn.diagram)
#attr(LABELS, "intersections") for groups of genes shown in venn


#4-16-19
##Sorted and indexed bam files so can pull specific sequence ranges
##Want to identify what fragments of what was identifed by alignment as MALAT1 (or Neat1, etc) are present (determines miRNA binding among other things)