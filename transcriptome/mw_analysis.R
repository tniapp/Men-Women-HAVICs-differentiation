## This code reproduces the analysis of RNA sequencing data from the article
## "Sex-Specific Gene Expression in Osteogenic Differentiation of Human Aortic Valve Interstitial Cells: Insights into Calcific Aortic Valve Disease"
## Code author: Polina Kuchur

library(readr)
library(DESeq2)
library(readxl)
library(org.Hs.eg.db)
library(vsn)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(EnhancedVolcano)
library(pheatmap)
library(DESeq2)
library(tidyverse)
library(ggvenn)
library(clusterProfiler)
library(VennDiagram)
library(ReactomePA)

setwd("./transcriptome")

##-- Compare male HAVICs on 10 days (differentiated) and 5 days (control) of differentiation

counts <- as.matrix(read.delim("./fc_M_CD.txt", header = TRUE, sep = "\t",row.names = 1))
head(counts)

# remove X before sample name
colnames(counts) <- substr(colnames(counts), 2, length(colnames(counts))) 
head(counts)

# upload sample's metadata
metadata <- read_xlsx("./MW_metadata.xlsx", sheet = 3)

colData <- data.frame(row.names = "sample",
                      sample = metadata$sample,
                      status = metadata$status,
                      gender = metadata$gender,
                      donor = metadata$donor)

# check the order of samples in colData and counts
all(rownames(colData) %in% colnames(counts))  # should be TRUE

colData <- colData[-c(11,12,13,14,15,16),]  # remove outlier samples

#get samples from common dataset based on their names
counts <- counts[,rownames(colData)]

# make sample features as factors
colData$status <- as.factor(colData$status)
colData$donor <- as.factor(colData$donor)

# prepare data for DESeq2 analysis
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = colData,
                                      design = ~donor+status)

# check and set the reference level for status (control)
dds$status <- relevel(dds$status, ref = "control")

# remove low-counted genes
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10,]
nrow(dds)

# remove suffix after dot from gene ensg name (it takes time)
for(i in rownames(dds)) {
  rownames(dds)[rownames(dds) == i] <- gsub(pattern = "\\.(.*)", "", i)
}
rownames(dds)

# add gene symbols to gene properties (ensg)
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)

## Start DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons

# extract DESeq2 results for status
res <- results(dds, contrast=c("status","differentiated","control"))

#summary of differentially expressed genes (degs)
summary(res, alpha = 0.05)

# MA-plot
plotMA(res)
dev.off()

# shrink log fold changes association with status
res_lfc <- lfcShrink(dds, coef="status_differentiated_vs_control", type="apeglm")

# MA-plot
plotMA(res_lfc)
dev.off()

#summary of differentially expressed genes (degs)
summary(res_lfc, alpha = 0.05)

# export DESeq2 results in the file
# write.csv(as.data.frame(res[order(res$padj),] ), file="./M_CD_degs.csv")

# Volcano plot for DEGs
# tiff(file="./M_CD_volcano.tiff",
#      units = "in",
#      width = 12,
#      height = 12,
#      res=300,
#      compression = "lzw")
EnhancedVolcano(res_lfc,
                lab = rownames(res_lfc),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-5, 6),
                ylim = c(0, 190),
                title ="Transcriptome of men VICs on 10 day vs 5 day of differentiation",
                labSize = 7.0,
                legendLabSize = 20.0,
                axisLabSize =20,
                titleLabSize = 20,
                subtitleLabSize = 10,
                captionLabSize = 20, 
                pointSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                col=c('grey','#CBD5E8','#B3E2CD','#FDCDAC'),
                selectLab = c("CLDN1", "COL1A1", "FN1", "ITGB3", "MMP2",
                              "MMP9", "ACE", "AXIN2", "CD34", "VCAM1",
                              "RNASE1", "FKBP5", "MAOA", "GLUL", "CARMIL1",
                              "IRS2", "METTL7A", "FBLN1", "FGD4", "DIAPH2",
                              "NEXN", "DNAJB4", "ALCAM", "SORBS1", "TEX2",
                              "NUDCD1", "IGFBP2", "SLC27A3", "GGT5",
                              "OSBPL10", "WDR59", "ADAMTSL4", "MAP2", 
                              "NOL3", "NID2", "TUT7", "STC2", "COL14A1",
                              "P4HA2", "CES2", "C4A", "CEMIP", "DAB2", 
                              "VASN", "ITIH3", "CCDC80", "GPC1", "EXT2",
                              "MMP14", "GPX4", "SERPINE2", "IGF2R", "ENG",
                              "MMP2", "IGFBP7", "ATP13A3", "IFNGR1", 
                              "SQSTM1", "GPX1", "PCDHGC3", "IFIT3", 
                              "CDH6", "CRIM1", "COL2A1", "BGN", "OSBPL8",
                              "TMEM59", "CA12", "LPIN1", "UBE2H", 
                              "ANGPTL2", "PAM", "PDGFRA", "GGA2", "NES",
                              "GBP2", "ACE", "ANPEP", "APCDD1", "CD86",
                              "CD36", "CD160", "CDH5", "F3", "KDR", 
                              "LYVE1", "SELP", "SOD3", "PTGIS"))
dev.off()

# normalize the data
rlt <- rlog(dds, blind=TRUE)

# inspect batch on PCA plot
pcaData <- plotPCA(rlt, intgroup=c("status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
pcaData[,5] <- c("c", "d", "c", "d",
                 "c", "d", "c", "d",
                 "c", "d",
                 "c", "d", "c", "d")
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")

# tiff(file="./M_CD_pca.tiff",
#      units = "in",
#      width = 9,
#      height = 7,
#      res=300,
#      compression = "lzw")
ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_colour_manual(values=brewer.pal(n = 9, "Paired")[c(2,8,4)])
dev.off()

## remove batch effect
assay(rlt) <- limma::removeBatchEffect(assay(rlt),
                                       batch = colData(dds)[,'donor'])

# Pheatmap for top DEGs
df <- as.data.frame(colData(dds)[,c("status")])
rownames(df) <- colnames(dds)
colnames(df) <- "cells"

res_lfc <- as.data.frame(res_lfc[order(res_lfc$padj),])
# select genes with padj < 0.05 
res_lfc <- subset(res_lfc, res_lfc$padj < 0.05)
# order by decreasing logfc
res_lfc <- res_lfc[order(abs(res_lfc$log2FoldChange), decreasing = TRUE),]

#select top genes
de<- rownames(res_lfc)[1:30]
de_mat <- assay(rlt)[de,]

# tiff(file="./M_CD_pheatmap.tiff",
#      units = "in",
#      width = 7,
#      height = 6,
#      res=300,
#      compression = "lzw")
pheatmap(t(scale(t(de_mat))),
         show_rownames = T,
         show_colnames = T,
         annotation_col = df,
         fontsize = 10)
dev.off()


## Enrichment analysis
vic_m_up <- subset(res_lfc, res_lfc$log2FoldChange > 1 & res_lfc$padj < 0.05 & !is.na(res_lfc$padj))
vic_m_up_genes <- rownames(vic_m_up)

vic_m_down <- subset(res_lfc, res_lfc$log2FoldChange < -1 & res_lfc$padj < 0.05 & !is.na(res_lfc$padj))
vic_m_down_genes <- rownames(vic_m_down)

# convert gene names into entrez ids
vic_m_up_genes.entrez <- clusterProfiler::bitr(vic_m_up_genes,
                                              fromType = "SYMBOL",
                                              toType = "ENTREZID", 
                                              OrgDb = org.Hs.eg.db)
vic_m_up_genes.entrez <- unique(vic_m_up_genes.entrez)

vic_m_down_genes.entrez <- clusterProfiler::bitr(vic_m_down_genes,
                                               fromType = "SYMBOL",
                                               toType = "ENTREZID", 
                                               OrgDb = org.Hs.eg.db)
vic_m_down_genes.entrez <- unique(vic_m_down_genes.entrez)

# create clusters
clusters <- list(men_dif_up = vic_m_up_genes.entrez$ENTREZID,
                 men_dif_down = vic_m_down_genes.entrez$ENTREZID)

## enrichGO
ck_go <- compareCluster(geneCluster = clusters, 
                     fun = enrichGO, 
                     OrgDb=org.Hs.eg.db, 
                     ont = "BP")
ck_pathways <- compareCluster(geneCluster = clusters, 
                              fun = enrichKEGG,
                              organism = "hsa")

# visualize these results
p1 <- dotplot(ck_go, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 15,
              title="DEGs in male HAVICs (GO: Biological processes)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

p2 <- dotplot(ck_pathways, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 15,
              title="DEGs in male HAVICs (KEGG pathways)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

# save the figure
# tiff(file="./M_CD_enrichment.tiff",
#      units = "in",
#      width = 13,
#      height = 9,
#      res=300)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()




##-- Compare female HAVICs on 10 days (differentiated) and 5 days (control) of differentiation

counts <- as.matrix(read.delim("./fc_W_CD.txt", header = TRUE, sep = "\t",row.names = 1))
head(counts)

# remove X before sample name
colnames(counts) <- substr(colnames(counts), 2, length(colnames(counts))) 
head(counts)

# sample info
metadata <- read_xlsx("./MW_metadata.xlsx", sheet = 4)

colData <- data.frame(row.names = "sample",
                      sample = metadata$sample,
                      status = metadata$status,
                      gender = metadata$gender,
                      donor = metadata$donor)

# check the order of samples in colData and counts
all(rownames(colData) %in% colnames(counts))  # should be TRUE

#get samples from common dataset based on their names
counts <- counts[,rownames(colData)]

# make sample features as factors
colData$status <- as.factor(colData$status)
colData$donor <- as.factor(colData$donor)

# prepare data for DESeq2 analysis
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = colData,
                                      design = ~donor+status)

# check and set the reference level for status (control)
dds$status <- relevel(dds$status, ref = "control")

# remove low-counted genes
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10,]
nrow(dds)

# remove suffix after dot from gene ensg name (it takes time)
for(i in rownames(dds)) {
  rownames(dds)[rownames(dds) == i] <- gsub(pattern = "\\.(.*)", "", i)
}
rownames(dds)

# add gene symbols to gene properties (ensg)
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)

## Start DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons

# extract DESeq2 results for status
res <- results(dds, contrast=c("status","differentiated","control"))

#summary of differentially expressed genes (degs)
summary(res, alpha = 0.05)

# MA-plot
plotMA(res)
dev.off()

# shrink log fold changes association with status
res_lfc <- lfcShrink(dds, coef="status_differentiated_vs_control", type="apeglm")

# MA-plot
plotMA(res_lfc)
dev.off()

#summary of differentially expressed genes (degs)
summary(res_lfc, alpha = 0.05)

# export DESeq2 results in the file
# write.csv(as.data.frame(res[order(res$padj),] ), file="./W_CD_degs.csv")

# Volcano plot
# tiff(file="./W_CD_volcano_scaled.tiff",
#      units = "in",
#      width = 12,
#      height = 12,
#      res=300,
#      compression = "lzw")
EnhancedVolcano(res_lfc,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-5, 10),
                ylim = c(0, 250),
                title ="Transcriptome of women VICs on 10 day vs 5 day of differentiation",
                labSize = 7.0,
                legendLabSize = 20.0,
                axisLabSize =20,
                titleLabSize = 20,
                subtitleLabSize = 10,
                captionLabSize = 20, 
                pointSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                col=c('grey','#CBD5E8','#B3E2CD','#FDCDAC'),
                selectLab = c("CLDN1", "COL1A1", "FN1", "ITGB3",
                              "MMP9", "FKBP5", "MAOA", "GLUL", "CARMIL1",
                              "IRS2", "METTL7A", "FBLN1", "FGD4", "DIAPH2",
                              "NEXN", "DNAJB4", "ALCAM", "SORBS1",
                              "NUDCD1", "IGFBP2", "SLC27A3", "GGT5",
                              "OSBPL10", "WDR59", "ADAMTSL4", "MAP2", 
                              "NOL3", "NID2", "TUT7", "STC2", "COL14A1",
                              "P4HA2", "CES2", "C4A", "CEMIP",  
                              "VASN", "ITIH3", "CCDC80", "GPC1", "EXT2",
                              "MMP14", "GPX4", "IGF2R", "ENG",
                              "IGFBP7", "ATP13A3", "IFNGR1", "SQSTM1", 
                              "GPX1", "PCDHGC3", "IFIT3", "CDH6", "CRIM1",
                              "COL2A1", "BGN", "OSBPL8", "TMEM59", "CA12",
                              "LPIN1", "UBE2H", "PAM", "PDGFRA",
                              "GGA2", "NES", "GBP2", "ACE", "ANPEP", 
                              "APCDD1", "AXIN2", "CD34", "CD36", "F3", 
                              "KDR", "LYVE1", "STAB2", "VCAM1", 
                              "RNASE1", "SOD3"))
dev.off()

# normalize the data
rlt <- rlog(dds, blind=TRUE)

## remove batch effect
assay(rlt) <- limma::removeBatchEffect(assay(rlt),
                                       batch = colData(dds)[,'donor'])

# inspect batch removal on PCA plot
# tiff(file="./W_CD_pca.tiff",
#      units = "in",
#      width = 9,
#      height = 7,
#      res=300,
#      compression = "lzw")
pcaData <- plotPCA(rlt, intgroup=c("status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
pcaData[,5] <- c("c", "d", "c", "d",
                 "c", "d", "c", "d",
                 "c", "d", "c", "d",
                 "c", "d", "c", "d",
                 "c", "d", "c", "d",
                 "c", "d")
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")
ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_colour_manual(values=brewer.pal(n = 9, "Paired")[c(2,8,4)])
dev.off()


# Pheatmap
df <- as.data.frame(colData(dds)[,c("status")])
rownames(df) <- colnames(dds)
colnames(df) <- "cells"

res_lfc <- as.data.frame(res_lfc[order(res_lfc$padj),])
# select genes with padj < 0.05 
res_lfc <- subset(res_lfc, res_lfc$padj < 0.05)
# order by decreasing logfc
res_lfc <- res_lfc[order(abs(res_lfc$log2FoldChange), decreasing = TRUE),]

#select top genes
de<- rownames(res_lfc)[1:30]
de_mat <- assay(rlt)[de,]

# tiff(file="./W_CD_pheatmap.tiff",
#      units = "in",
#      width = 7,
#      height = 6,
#      res=300,
#      compression = "lzw")
pheatmap(t(scale(t(de_mat))),
         show_rownames = T,
         show_colnames = T,
         annotation_col = df,
         fontsize = 10)
dev.off()


## Enrichment analysis
vic_w_up <- subset(res_lfc, res_lfc$log2FoldChange > 1 & res_lfc$padj < 0.05 & !is.na(res_lfc$padj))
vic_w_up_genes <- rownames(vic_w_up)

vic_w_down <- subset(res_lfc, res_lfc$log2FoldChange < -1 & res_lfc$padj < 0.05 & !is.na(res_lfc$padj))
vic_w_down_genes <- rownames(vic_w_down)


# convert gene names into entrez ids
vic_w_up_genes.entrez <- clusterProfiler::bitr(vic_w_up_genes,
                                               fromType = "SYMBOL",
                                               toType = "ENTREZID", 
                                               OrgDb = org.Hs.eg.db)
vic_w_up_genes.entrez <- unique(vic_w_up_genes.entrez)

vic_w_down_genes.entrez <- clusterProfiler::bitr(vic_w_down_genes,
                                                 fromType = "SYMBOL",
                                                 toType = "ENTREZID", 
                                                 OrgDb = org.Hs.eg.db)
vic_w_down_genes.entrez <- unique(vic_w_down_genes.entrez)


# create clusters
clusters <- list(women_dif_up = vic_w_up_genes.entrez$ENTREZID,
                 women_dif_down = vic_w_down_genes.entrez$ENTREZID)

## enrichGO
ck_go <- compareCluster(geneCluster = clusters, 
                        fun = enrichGO, 
                        OrgDb=org.Hs.eg.db, 
                        ont = "BP")
ck_pathways <- compareCluster(geneCluster = clusters, 
                              fun = enrichKEGG,
                              organism = "hsa")


# visualize these results
p1 <- dotplot(ck_go, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="DEGs in female VICs (Biological processes)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

p2 <- dotplot(ck_pathways, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="DEGs in female VICs (KEGG pathways)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

# save the figure
# tiff(file="./W_CD_enrichment.tiff",
     units = "in",
     width = 13,
     height = 8,
     res=300)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()



### Check for overlap between men and woman upregulated DEGs
m <- read.csv("./M_CD_degs.csv", row.names = 1)
m_deg <- subset(m, m$log2FoldChange > 1 & m$padj < 0.05 & !is.na(m$padj))
m_genes <- rownames(m_deg)

w <- read.csv("./W_CD_degs.csv", row.names = 1)
w_deg <- subset(w, w$log2FoldChange > 1 & w$padj < 0.05 & !is.na(w$padj))
w_genes <- rownames(w_deg)

venn.diagram(
  x = list(m_genes,
           w_genes),
  category.names = c("M 10d dif up", 
                     "W 10d dif up"),
  filename = './M&W_10d_up_DEG_venn.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 700 , 
  width = 700 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 1,
  lty = 1,
  fill = c("#377EB8", "#FF7F00"),
  main.cex = 10,
  
  # Numbers
  cex = 1,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.pos = c(200, 10),
  cat.dist = c(0.05,0.05),
  cat.fontfamily = "sans"
)


w_unique <- setdiff(w_genes, m_genes)  # up DEGs unique for woman HAVIC's on 10 days of differentiation
m_unique <- setdiff(m_genes, w_genes)  # up DEGs unique for men HAVIC's on 10 days of differentiation


## enrichment of unique genes
# convert gene names into entrez ids
w_unique.entrez <- clusterProfiler::bitr(w_unique,
                                         fromType = "SYMBOL",
                                         toType = "ENTREZID",
                                         OrgDb = org.Hs.eg.db)
w_unique.entrez <- unique(w_unique.entrez)

m_unique.entrez <- clusterProfiler::bitr(m_unique,
                                         fromType = "SYMBOL",
                                         toType = "ENTREZID",
                                         OrgDb = org.Hs.eg.db)
m_unique.entrez <- unique(m_unique.entrez)


# create clusters
clusters <- list(women_dif_up_unique = w_unique.entrez$ENTREZID,
                 men_dif_up_unique = m_unique.entrez$ENTREZID)

## enrichGO
ck_go <- compareCluster(geneCluster = clusters, 
                        fun = enrichGO, 
                        OrgDb=org.Hs.eg.db, 
                        ont = "BP")

# visualize these results
p1 <- dotplot(ck_go, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="up DEGs unique for W 10d dif or M 10d dif (GO: Biological processes)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

# save the figure
# tiff(file="./W_M_up_10d_unique_enrichment.tiff",
#      units = "in",
#      width = 8,
#      height = 7,
#      res=300)
gridExtra::grid.arrange(p1, ncol = 1)
dev.off()



##-- Compare women HAVICs (10 days of diff) and men HAVICs (10 days of diff)
counts <- as.matrix(read.delim("./fc_D_MW.txt", header = TRUE, sep = "\t",row.names = 1))
head(counts)

# remove X before sample name
colnames(counts) <- substr(colnames(counts), 2, length(colnames(counts))) 
head(counts)

# sample info
metadata <- read_xlsx("./MW_metadata.xlsx", sheet = 2)

colData <- data.frame(row.names = "sample",
                      sample = metadata$sample,
                      status = metadata$status,
                      gender = metadata$gender,
                      donor = metadata$donor)

# check the order of samples in colData and counts
all(rownames(colData) %in% colnames(counts))  # should be TRUE

colData <- colData[-c(13,17,18),]

#get samples from common dataset based on their names
counts <- counts[,rownames(colData)]

# make sample features as factors
colData$gender <- as.factor(colData$gender)


# prepare data for DESeq2 analysis
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = colData,
                                      design = ~gender)

# check and set the reference level for status (woman)
dds$gender <- relevel(dds$gender, ref = "woman")

# remove low-counted genes
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10,]
nrow(dds)


# remove suffix after dot from gene ensg name (it takes time)
for(i in rownames(dds)) {
  rownames(dds)[rownames(dds) == i] <- gsub(pattern = "\\.(.*)", "", i)
}
rownames(dds)

# add gene symbols to gene properties (ensg)
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)

## Start DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons

# extract DESeq2 results for status
res <- results(dds, contrast=c("gender","man","woman"))

#summary of differentially expressed genes (degs)
summary(res, alpha = 0.05)

# MA-plot
plotMA(res)
dev.off()

# shrink log fold changes association with status
res_lfc <- lfcShrink(dds, coef="gender_man_vs_woman", type="apeglm")

# MA-plot
plotMA(res_lfc)
dev.off()

# summary of differentially expressed genes (degs)
summary(res_lfc, alpha = 0.05)

# export DESeq2 results in the file
# write.csv(as.data.frame(res[order(res$padj),] ), file="./W_D_vs_M_D_degs.csv")

# Volcano plot
# tiff(file="./WM_D_volcano.tiff",
#      units = "in",
#      width = 9,
#      height = 7,
#      res=300,
#      compression = "lzw")
EnhancedVolcano(res_lfc,
                lab = rownames(res_lfc),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title ="VIC male 10d dif vs female 10d dif",
                labSize = 3,
                boxedLabels = F,
                col=c('black',
                             '#CBD5E8', 
                             '#B3E2CD', 
                             '#FDCDAC'),
                             colAlpha = 1)
dev.off()


# normalize the data
rlt <- rlog(dds, blind=TRUE)

# inspect batch removal on PCA plot
# tiff(file="./WM_D_pca.tiff",
#      units = "in",
#      width = 9,
#      height = 7,
#      res=300,
#      compression = "lzw")
pcaData <- plotPCA(rlt, intgroup=c("gender"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
pcaData[,5] <- c("d", "d", "d", "d", "d",
                 "d", "d", "d", "d", "d",
                 "d", "d", "d", "d", "d",
                 "d", "d", "d")
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")
ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_colour_manual(values=brewer.pal(n = 9, "Paired")[c(2,8,4)])
dev.off()


# Pheatmap
df <- as.data.frame(colData(dds)[,c("gender")])
rownames(df) <- colnames(dds)
colnames(df) <- "cells"

res_lfc <- as.data.frame(res_lfc[order(res_lfc$padj),])
# select genes with padj < 0.05 
res_lfc <- subset(res_lfc, res_lfc$padj < 0.05 & abs(res_lfc$log2FoldChange) > 1)
# order by decreasing logfc
res_lfc <- res_lfc[order(abs(res_lfc$log2FoldChange), decreasing = TRUE),]

#select top genes
de<- rownames(res_lfc)[1:40]
de_mat <- assay(rlt)[de,]

# tiff(file="./WM_D_pheatmap.tiff",
#      units = "in",
#      width = 7,
#      height = 6,
#      res=300,
#      compression = "lzw")
pheatmap(t(scale(t(de_mat))),
         show_rownames = T,
         show_colnames = T,
         annotation_col = df,
         fontsize = 10)
dev.off()


## Enrichment analysis
vic_m_up <- subset(res_lfc, res_lfc$log2FoldChange > 1 & res_lfc$padj < 0.05 & !is.na(res_lfc$padj))
vic_m_up_genes <- rownames(vic_m_up)

vic_w_up <- subset(res_lfc, res_lfc$log2FoldChange < -1 & res_lfc$padj < 0.05 & !is.na(res_lfc$padj))
vic_w_up_genes <- rownames(vic_w_up)


# convert gene names into entrez ids
vic_m_up_genes.entrez <- clusterProfiler::bitr(vic_m_up_genes,
                                               fromType = "SYMBOL",
                                               toType = "ENTREZID", 
                                               OrgDb = org.Hs.eg.db)
vic_m_up_genes.entrez <- unique(vic_m_up_genes.entrez)

vic_w_up_genes.entrez <- clusterProfiler::bitr(vic_w_up_genes,
                                                 fromType = "SYMBOL",
                                                 toType = "ENTREZID", 
                                                 OrgDb = org.Hs.eg.db)
vic_w_up_genes.entrez <- unique(vic_w_up_genes.entrez)


# create clusters
clusters <- list(men_dif_up = vic_m_up_genes.entrez$ENTREZID,
                 women_dif_up = vic_w_up_genes.entrez$ENTREZID)

## enrichGO
ck_go <- compareCluster(geneCluster = clusters, 
                        fun = enrichGO, 
                        OrgDb=org.Hs.eg.db, 
                        ont = "BP")
ck_pathways <- compareCluster(geneCluster = clusters, 
                              fun = enrichPathway,
                              organism = "human")

# visualize these results
p1 <- dotplot(ck_go, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="DEGs in 10d dif male HAVICs vs 10d dif female HAVICs (Biological processes)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

p2 <- dotplot(ck_pathways, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="DEGs in 10d dif male HAVICs vs 10d dif female HAVICs (Reactome pathways)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

# save the figure
# tiff(file="./WM_D_enrichment.tiff",
#      units = "in",
#      width = 13,
#      height = 8,
#      res=300)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()



##-- Compare female 5 days of diff HAVICs and male 5 days of diff HAVICs
counts <- as.matrix(read.delim("./fc_C_MW.txt", header = TRUE, sep = "\t",row.names = 1))
head(counts)

# remove X before sample name
colnames(counts) <- substr(colnames(counts), 2, length(colnames(counts))) 
head(counts)

# sample info
metadata <- read_xlsx("./MW_metadata.xlsx", sheet = 1)

colData <- data.frame(row.names = "sample",
                      sample = metadata$sample,
                      status = metadata$status,
                      gender = metadata$gender,
                      donor = metadata$donor)

# check the order of samples in colData and counts
all(rownames(colData) %in% colnames(counts))  # should be TRUE

colData <- colData[-c(13,17,18),] # remove three man donors

#get samples from common dataset based on their names
counts <- counts[,rownames(colData)]

# make sample features as factors
colData$gender <- as.factor(colData$gender)

# prepare data for DESeq2 analysis
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = colData,
                                      design = ~gender)

# check and set the reference level for status (woman)
dds$gender <- relevel(dds$gender, ref = "woman")

# remove low-counted genes
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10,]
nrow(dds)

# remove suffix after dot from gene ensg name (it takes time)
for(i in rownames(dds)) {
  rownames(dds)[rownames(dds) == i] <- gsub(pattern = "\\.(.*)", "", i)
}
rownames(dds)


# add gene symbols to gene properties (ensg)
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}

# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)

## Start DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons

# extract DESeq2 results for status
res <- results(dds, contrast=c("gender","man","woman"))

#summary of differentially expressed genes (degs)
summary(res, alpha = 0.05)

# MA-plot
plotMA(res)
dev.off()

# shrink log fold changes association with status
res_lfc <- lfcShrink(dds, coef="gender_man_vs_woman", type="apeglm")

# MA-plot
plotMA(res_lfc)
dev.off()

#summary of differentially expressed genes (degs)
summary(res_lfc, alpha = 0.05)

# export DESeq2 results in the file
# write.csv(as.data.frame(res[order(res$padj),] ), file="./MW_C_degs.csv")

# Volcano plot
# tiff(file="./MW_C_volcano.tiff",
#      units = "in",
#      width = 9,
#      height = 7,
#      res=300,
#      compression = "lzw")
EnhancedVolcano(res_lfc,
                lab = rownames(res_lfc),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                title ="HAVICs male (5 days) vs female (5 days)",
                labSize = 3,
                boxedLabels = F,
                col=c('black',
                             '#CBD5E8', 
                             '#B3E2CD', 
                             '#FDCDAC'),
                             colAlpha = 1)
dev.off()


# normalize the data
rlt <- rlog(dds, blind=TRUE)

# inspect batch removal on PCA plot
# tiff(file="./MW_C_pca.tiff",
#      units = "in",
#      width = 9,
#      height = 7,
#      res=300,
#      compression = "lzw")
pcaData <- plotPCA(rlt, intgroup=c("gender"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
pcaData[,5] <- c("c", "c", "c", "c", "c",
                 "c", "c", "c", "c", "c",
                 "c", "c", "c", "c", "c",
                 "c", "c", "c")
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")
ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_colour_manual(values=brewer.pal(n = 9, "Paired")[c(2,8,4)])
dev.off()


# Pheatmap
df <- as.data.frame(colData(dds)[,c("gender")])
rownames(df) <- colnames(dds)
colnames(df) <- "cells"

res_lfc <- as.data.frame(res_lfc[order(res_lfc$padj),])
# select genes with padj < 0.05 
res_lfc <- subset(res_lfc, res_lfc$padj < 0.05 & abs(res_lfc$log2FoldChange) > 1)
# order by decreasing logfc
res_lfc <- res_lfc[order(abs(res_lfc$log2FoldChange), decreasing = TRUE),]

#select top genes
de<- rownames(res_lfc)
de_mat <- assay(rlt)[de,]

# tiff(file="./MW_C_pheatmap.tiff",
#      units = "in",
#      width = 7,
#      height = 6,
#      res=300,
#      compression = "lzw")
pheatmap(t(scale(t(de_mat))),
         show_rownames = T,
         show_colnames = T,
         annotation_col = df,
         fontsize = 10)
dev.off()


## Enrichment analysis
vic_m_up <- subset(res_lfc, res_lfc$log2FoldChange > 1 & res_lfc$padj < 0.05 & !is.na(res_lfc$padj))
vic_m_up_genes <- rownames(vic_m_up)

vic_w_up <- subset(res_lfc, res_lfc$log2FoldChange < -1 & res_lfc$padj < 0.05 & !is.na(res_lfc$padj))
vic_w_up_genes <- rownames(vic_w_up)


# convert gene names into entrez ids
vic_m_up_genes.entrez <- clusterProfiler::bitr(vic_m_up_genes,
                                               fromType = "SYMBOL",
                                               toType = "ENTREZID", 
                                               OrgDb = org.Hs.eg.db)
vic_m_up_genes.entrez <- unique(vic_m_up_genes.entrez)

vic_w_up_genes.entrez <- clusterProfiler::bitr(vic_w_up_genes,
                                               fromType = "SYMBOL",
                                               toType = "ENTREZID", 
                                               OrgDb = org.Hs.eg.db)
vic_w_up_genes.entrez <- unique(vic_w_up_genes.entrez)


# create clusters
clusters <- list(men_dif_up = vic_m_up_genes.entrez$ENTREZID,
                 women_dif_up = vic_w_up_genes.entrez$ENTREZID)

## enrichGO
ck_go <- compareCluster(geneCluster = clusters, 
                        fun = enrichGO, 
                        OrgDb=org.Hs.eg.db, 
                        ont = "BP")
ck_pathways <- compareCluster(geneCluster = clusters, 
                              fun = enrichPathway,
                              organism = "human")

# visualize these results
p1 <- dotplot(ck_go, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="DEGs for male HAVICs (5 days) vs female HAVICs (5 days) (Biological processes)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

p2 <- dotplot(ck_pathways, 
              showCategory = 10,
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="DEGs in contr male VICs vs contr female VICs (Reactome pathways)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12))

# save the figure
# tiff(file="./MW_C_enrichment.tiff",
#      units = "in",
#      width = 13,
#      height = 8,
#      res=300)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()




#---------------------- WNT --------------------------#
setwd("../")

counts <- as.matrix(read.delim("./wnt/wnt_counts.csv", header = TRUE, sep = ";",row.names = 1))
head(counts)

# sample info
metadata <- read_xlsx("./wnt/metadata_wnt.xlsx", sheet = 1)

colData <- data.frame(row.names = "sample",
                      sample = metadata$sample,
                      status = metadata$status,
                      replicate = metadata$replicate)

# check the order of samples in colData and counts
all(rownames(colData) %in% colnames(counts))  # should be TRUE

#get samples from common dataset based on their names
counts <- counts[,rownames(colData)]

# make sample features as factors
colData$status <- as.factor(colData$status)
colData$replicate <- as.factor(colData$replicate)

# prepare data for DESeq2 analysis
dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                      colData = colData,
                                      design = ~status)

dds <- collapseReplicates(dds, dds$replicate)

# check and set the reference level for status (control)
dds$status <- relevel(dds$status, ref = "empty")

# remove low-counted genes
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 10,]
nrow(dds)

# remove suffix after dot from gene ensg name (it takes time)
for(i in rownames(dds)) {
  rownames(dds)[rownames(dds) == i] <- gsub(pattern = "\\.(.*)", "", i)
}
rownames(dds)

# add gene symbols to gene properties
genes <- rownames(dds)
symbols <- mapIds(org.Hs.eg.db, keys = genes,
                  column = c('SYMBOL'),
                  keytype = 'ENSEMBL')
genes_symbol <- data.frame(genes, symbols)

counter=1
for(gene in genes_symbol$symbols){
  if (is.na(gene)) {
    genes_symbol[counter,2] <- genes_symbol[counter,1]
  }
  counter = counter+1
}
# change ensg to gene symbols where possible
rownames(dds) <- genes_symbol$symbols
nrow(dds)

## Start DESeq2 analysis
dds <- DESeq(dds)
resultsNames(dds)  # see all comparisons

res <- results(dds, contrast=c("status","Wnt","empty"))
summary(res)

# export DESeq results in the file
# write.csv(as.data.frame(res[order(res$padj),] ), file="./wnt/wnt_vs_no_empty_deseq.csv")

# MA-plot
plotMA(res)
dev.off()

# shrink log fold changes association with status
res_lfc <- lfcShrink(dds, coef="status_Wnt_vs_empty", type="apeglm")

# MA-plot
plotMA(res_lfc)
dev.off()

#summary of differentially expressed genes (degs)
summary(res_lfc, alpha = 0.05)
res <- res_lfc

top_degs <- as.data.frame(subset(res, abs(res$log2FoldChange) > 1 & res$padj < 0.05 & !is.na(res$padj))) %>% arrange(desc(padj))
top20_degs <- tail(rownames(top_degs), n = 20)

# tiff(file="./wnt/wnt_volcano.tiff",
#      units = "in",
#      width = 15,
#      height = 10,
#      res=300,
#      compression = "lzw")
EnhancedVolcano(res_lfc,
                lab = rownames(res_lfc),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-12, 15),
                title ="Transcriptome of VICs with activated Wnt",
                labSize = 7.0,
                legendLabSize = 20.0,
                axisLabSize =20,
                titleLabSize = 20,
                subtitleLabSize = 10,
                captionLabSize = 20, 
                pointSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = F,
                colAlpha = 4/5,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                col=c('grey','#CBD5E8','#B3E2CD','#FDCDAC'),
                selectLab = c(top20_degs, "APCDD1",
                              "AXIN2", "PRKCH", "RNASE1",
                              "BMP6", "TIMP3", "MGP",
                              "ACKR3"))
dev.off()


# normalize the data
rlt <- rlog(dds, blind=TRUE)

# inspect batch on PCA plot
pcaData <- plotPCA(rlt, intgroup=c("status"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
names <- pcaData[,5]
colnames(pcaData) <-c("PC1","PC2","group","condition","replicate")

# tiff(file="./wnt/wnt_pca.tiff",
#      units = "in",
#      width = 5,
#      height = 5,
#      res=300,
#      compression = "lzw")
ggplot(pcaData, aes(PC1, PC2, shape=replicate, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  scale_colour_manual(values=brewer.pal(n = 9, "Paired")[c(2,8,4)])
dev.off()


# Pheatmap
df <- as.data.frame(colData(dds)[,c("status")])
rownames(df) <- colnames(dds)
colnames(df) <- "cells"

res_lfc <- as.data.frame(res_lfc[order(res_lfc$padj),])
# select genes with padj < 0.05 
res_lfc <- subset(res_lfc, res_lfc$padj < 0.05)

#select top genes
de<- rownames(res_lfc)[1:30]
de_mat <- assay(rlt)[de,]

# tiff(file="./wnt/wnt_heatmap.tiff",
#      units = "in",
#      width = 7,
#      height = 8,
#      res=300,
#      compression = "lzw")
pheatmap(t(scale(t(de_mat))),
         show_rownames = T,
         show_colnames = T,
         annotation_col = df,
         fontsize = 14)
dev.off()


## Enrichment analysis
deg <- read.csv("./wnt/wnt_vs_no_empty_deseq.csv")

up <- subset(deg, deg$log2FoldChange > 1 & deg$padj < 0.05 & !is.na(deg$padj))
up_genes <- up$X

down <- subset(deg, deg$log2FoldChange < -1 & deg$padj < 0.05 & !is.na(deg$padj))
down_genes <- down$X


# convert gene names into entrez ids
up_genes.entrez <- clusterProfiler::bitr(up_genes,
                                               fromType = "SYMBOL",
                                               toType = "ENTREZID", 
                                               OrgDb = org.Hs.eg.db)
up_genes.entrez <- unique(up_genes.entrez)

down_genes.entrez <- clusterProfiler::bitr(down_genes,
                                         fromType = "SYMBOL",
                                         toType = "ENTREZID", 
                                         OrgDb = org.Hs.eg.db)
down_genes.entrez <- unique(down_genes.entrez)


# create clusters
clusters <- list("WNT↑" = up_genes.entrez$ENTREZID,
                 "WNT↓" = down_genes.entrez$ENTREZID)

## enrichGO
ck_go <- compareCluster(geneCluster = clusters, 
                        fun = enrichGO,
                        minGSSize = 5,
                        OrgDb=org.Hs.eg.db, 
                        ont = "BP")
ck_pathways <- compareCluster(geneCluster = clusters, 
                              fun = enrichKEGG,
                              minGSSize = 2,
                              organism = "hsa")


# visualize these results
p1 <- dotplot(ck_go,
              showCategory = c("cytoplasmic pattern recognition receptor signaling pathway", "positive regulation of interferon-beta production", "interferon-beta production", "	
positive regulation of cytokine production", "interleukin-27-mediated signaling pathway", "regulation of prostaglandin secretion", "positive regulation of prostaglandin secretion", "extracellular matrix organization", "extracellular structure organization", "	
response to calcium ion", "connective tissue development", "negative regulation of chondrocyte differentiation", "somatostatin receptor signaling pathway", "somatostatin signaling pathway", "negative regulation of protein kinase activity", "regulation of membrane potential", "negative regulation of cartilage development", "negative regulation of protein phosphorylation"),
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="WNT GO (Biological processes)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12), plot.title = element_text(hjust = 0.5))

p2 <- dotplot(ck_pathways, 
              showCategory = 30,
              by = "Count",
              color = "p.adjust",
              font.size = 12.5,
              title="WNT pathways (KEGG)",
              includeAll=TRUE) + scale_colour_distiller(palette="Set2") +
  theme(axis.text.y=element_text(size=12), plot.title = element_text(hjust = 0.5))

up_wiki <- enrichWP(gene = up_genes.entrez$ENTREZID,
                                organism = "Homo sapiens",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                                minGSSize = 5,
                                maxGSSize = 500,
                                qvalueCutoff = 0.2)

# tiff(file="./wnt/wnt_upregulated_enrichWiki.tiff",
#      units = "in",
#      width = 8,
#      height = 7,
#      res=300)
dotplot(up_wiki, showCategory = 30, font.size = 12, label_format = 50) + scale_colour_distiller(palette="Set2")
dev.off()

# save the figure
# tiff(file="./wnt/wnt_enrichBP_enrichKEGG.tiff",
#      units = "in",
#      width = 11,
#      height = 7.5,
#      res=300)
gridExtra::grid.arrange(p1, p2, ncol = 2)
dev.off()

