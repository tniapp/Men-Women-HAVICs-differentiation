## This code reproduces the analysis of proteomic data from the article
## "Sex-Specific Gene Expression in Osteogenic Differentiation of Human Aortic Valve Interstitial Cells: Insights into Calcific Aortic Valve Disease"
## Code author: Polina Kuchur

library(DEP)
library(readxl)
library(dplyr)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(mixOmics)

setwd("./proteome")

##------------------------------------------------------------------------##
########----- Compare men 10 days HAVICs with men 5 days HAVICs -----#######
##------------------------------------------------------------------------##

### upload protein counts before NA imputation
data <- read.csv("./Area_data_42_samples_full_spec_lib.csv", sep = ",")

#metadata
experimental_design <- read_xlsx("./Sample_data_days.xlsx", sheet = 2)

#subset samples of interest (men 5 days and 10 days)
data_men <- data[,c("Protein.Ids", experimental_design$label)]

data[data == 0] = NA

#filtering
data_men$na.count <- apply(data_men, 1, function(x) sum(is.na(x)))
data_men$na.count
nrow(data_men)

data_men <- subset(data_men, data_men$na.count < 5)
nrow(data_men)

#convert uniprot ids into gene names
proteins <- data_men$Protein.Ids 
symbols <- mapIds(org.Hs.eg.db, keys = proteins,
                  column = c('SYMBOL'),
                  keytype = 'UNIPROT')
proteins_symbol <- data.frame(proteins, symbols)

counter=1
for(protein in proteins_symbol$symbols){
  if (is.na(protein)) {
    proteins_symbol[counter,2] <- proteins_symbol[counter,1]
  }
  counter = counter+1
}

# change uniprot ids to protein symbols where possible
data_men$Protein <- proteins_symbol$symbols

#remove contaminants
data_men <- data_men[!grepl("contam_sp", data_men$Protein.Ids),]

#checking for duplicated proteins
data_men$Protein.Ids %>% duplicated() %>% any()

#prepare data for DEP analysis
data_unique <- make_unique(data_men, "Protein.Ids", "Protein", delim = ";")

#checking for duplicated proteins
data_unique$name %>% duplicated() %>% any()

#generate a SummarizedExperiment object using an experimental design
#and perform log2-transformation
ALM_columns <- grep("ALM", colnames(data_unique))
data_se <- make_se(data_unique, ALM_columns, experimental_design)

#plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)


#checking for missing values
data_filt <- filter_missval(data_se, thr = 1)
plot_frequency(data_filt)

#protein number
plot_numbers(data_se)
dev.off()

#normalize the data (if necessary)
data_norm <- normalize_vsn(data_filt)

#visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm) 

#plot a heatmap of proteins with missing values
plot_missval(data_filt)

#impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "knn")

#plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)


## Differential enrichment analysis  based on linear models and empherical Bayes statistics
#test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "X5d", design_formula = formula(~ 0 + condition))  # specify the MSCs type for the "control" option

#denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05)

#save dep results as a file
dep_results <- get_results(dep)
colnames(dep_results)
dep_results %>% filter(significant) %>% nrow()

#write.csv(x = dep_results, file = "M10d_vs_M5d_prot_dep.csv") 

dep$replicate <- ""

## PCA
# tiff(file="./M10d_VS_M5d_PCA.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4, label = T) + scale_color_brewer(palette = "Set2")
dev.off()


#PLS-DA 
ordination.optimum.splsda <- splsda(t(data_norm@assays@data@listData[[1]]), experimental_design$condition, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

# tiff('PLSDA_M10day_vs_M5d_comp.tiff', 
#      units="in", 
#      width=10, 
#      height=8, 
#      res=300, 
#      compression = 'lzw')
layout(matrix(c(1, 2)))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
dev.off()

# tiff('PLSDA_M10day_vs_M5d.tiff', 
#      units="in", 
#      width=10, 
#      height=8, 
#      res=300, 
#      compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 2, size.axis = 1.0, size.xlabel = 1.2, size.ylabel = 1.2, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()

rownames(dep) <- dep@elementMetadata@listData[["Protein"]]

## Heatmap of all significant proteins
# tiff(file="./M10d_VS_M5d_heatmap.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 5,
#     res=300)
plot_heatmap(dep, type = "centered",
             kmeans = TRUE,
             k = 2,
             col_limit = 4, 
             show_row_names = TRUE,
             indicate = c("condition"),
)
dev.off()

## Volcano plot for the contrast
# tiff(file="./M10d_VS_M5d_volcano_DEP.tiff",  
#     units = "in",
#     width = 12,
#     height = 12,
#     res=300)
data <- dep_results
data_p <- data[, c(2, 7, 3, 4)]
names(data_p)[1] <- "GeneName"
names(data_p)[2] <- "logFC"
names(data_p)[3] <- "P.Value"
names(data_p)[4] <- "adj.P.Val"

EnhancedVolcano(data_p,
                lab = data_p$GeneName,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-3, 3),
                ylim = c(0, 13),
                title ="Proteome of men VICs on 10 day vs 5 day of differentiation",
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
                selectLab = c("METTL7A", "ALCAM", "CARMIL1", 
                              "FKBP5", "MAOA", 
                              "FGD4", "IRS2", "GLUL",
                              "SORBS1", "COL14A1", 
                              "COL1A1", "MMP2","STAB2", "FBLN1", 
                              "DNAJB4", "TEX2", 
                              "P4HA2", "CES2", "C4A", "CEMIP", "DAB2", "VASN",
                              "ITIH3", "CCDC80", "GPC1", "EXT2", "MMP14", 
                              "GPX4", "SERPINE2", "ENG", "IGFBP7", 
                              "IFNGR1","GGA2"))
dev.off()



##----------------------------------------------------------##
########----- Compare women 10 days with women 5 days -----#######
##----------------------------------------------------------##

### upload protein counts before NA imputation
data <- read.csv("./Area_data_42_samples_full_spec_lib.csv", sep = ",")

#metadata
experimental_design <- read_xlsx("./Sample_data_days.xlsx", sheet = 3)

#subset samples of interest (women 5 days and 10 days)
data_women <- data[,c("Protein.Ids", experimental_design$label)]

data_women[data_women == 0] = NA

#filtering
data_women$na.count <- apply(data_women, 1, function(x) sum(is.na(x)))
data_women$na.count
nrow(data_women)

data_women <- subset(data_women, data_women$na.count < 5)
nrow(data_women)


#convert uniprot ids into gene names
proteins <- data_women$Protein.Ids 
symbols <- mapIds(org.Hs.eg.db, keys = proteins,
                  column = c('SYMBOL'),
                  keytype = 'UNIPROT')
proteins_symbol <- data.frame(proteins, symbols)

counter=1
for(protein in proteins_symbol$symbols){
  if (is.na(protein)) {
    proteins_symbol[counter,2] <- proteins_symbol[counter,1]
  }
  counter = counter+1
}

# change uniprot ids to protein symbols where possible
data_women$Protein <- proteins_symbol$symbols

#remove contaminants
data_women <- data_women[!grepl("contam_sp", data_women$Protein.Ids),]

#checking for duplicated proteins
data_women$Protein.Ids %>% duplicated() %>% any()

#prepare data for DEP analysis
data_women_unique <- make_unique(data_women, "Protein.Ids", "Protein", delim = ";")

#checking for duplicated proteins
data_women_unique$name %>% duplicated() %>% any()

#generate a SummarizedExperiment object using an experimental design
#and perform log2-transformation
ALM_columns <- grep("ALM", colnames(data_women_unique))
data_se <- make_se(data_women_unique, ALM_columns, experimental_design)

#plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)


#checking for missing values
data_filt <- filter_missval(data_se, thr = 1)
plot_frequency(data_filt)

#protein number
plot_numbers(data_se)
dev.off()

#normalize the data (if necessary)
data_norm <- normalize_vsn(data_filt)

#visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm) 

#plot a heatmap of proteins with missing values
plot_missval(data_filt)

#impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "knn")

#plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)


## Differential enrichment analysis  based on linear models and empherical Bayes statistics
#test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "X5d", design_formula = formula(~ 0 + condition))  # specify the MSCs type for the "control" option

#denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05)

#save dep results as a file
dep_results <- get_results(dep)
colnames(dep_results)
dep_results %>% filter(significant) %>% nrow()

#write.csv(x = dep_results, file = "W10d_vs_W5d_prot_dep.csv") 

dep$replicate <- ""

## PCA
# tiff(file="./W10d_VS_W5d_PCA.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4, label = T) + scale_color_brewer(palette = "Set2")
dev.off()


#PLS-DA 
ordination.optimum.splsda <- splsda(t(data_norm@assays@data@listData[[1]]), experimental_design$condition, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

# tiff('PLSDA_W10day_vs_W5d_comp.tiff',
#      units="in",
#      width=10,
#      height=8,
#      res=300,
#      compression = 'lzw')
layout(matrix(c(1, 2)))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
dev.off()

# tiff('PLSDA_W10day_vs_W5d.tiff',
#      units="in",
#      width=10,
#      height=8,
#      res=300,
#      compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 2, size.axis = 1.0, size.xlabel = 1.2, size.ylabel = 1.2, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()

rownames(dep) <- dep@elementMetadata@listData[["Protein"]]

## Heatmap of all significant proteins
# tiff(file="./W10d_VS_W5d_heatmap.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 7,
#     res=300)
plot_heatmap(dep, type = "centered",
             kmeans = TRUE,
             k = 2,
             col_limit = 4, 
             show_row_names = TRUE,
             indicate = c("condition"),
)
dev.off()

## Volcano plot for the contrast
# tiff(file="./W10d_VS_W5d_volcano_DEP.tiff",
#     units = "in",
#     width = 12,
#     height = 12,
#     res=300)
data <- dep_results
data_p <- data[, c(2, 7, 3, 4)]
names(data_p)[1] <- "GeneName"
names(data_p)[2] <- "logFC"
names(data_p)[3] <- "P.Value"
names(data_p)[4] <- "adj.P.Val"

EnhancedVolcano(data_p,
                lab = data_p$GeneName,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                xlim = c(-3, 3),
                ylim = c(0, 13),
                title ="Proteome of women VICs on 10 day vs 5 day of differentiation",
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
                selectLab = c("FKBP5", "METTL7A", "MAOA", 
                              "GLUL", "DIAPH2", "CCDC80", 
                              "ADAMTSL4", "ALCAM", "CARMIL1", 
                              "IRS2", "NEXN", "COL14A1", "FGD4", 
                              "GGT5", "DNAJB4", "SORBS1", "C4A",
                              "DKK1", "COL1A1", "MMP2", "IRS2", "METTL7A", 
                              "FGD4", "DIAPH2", "DNAJB4",
                              "GGT5", "PDGFRA",
                              "OSBPL10", "ADAMTSL4",
                              "P4HA2", "C4A", "CEMIP", "DAB2",
                              "VASN", "ITIH3", "CCDC80",
                              "MMP14", "GPX4", "SERPINE2", "ENG",
                              "IGFBP7", "ATP13A3", 
                              "PCDHGC3", "CRIM1", "COL2A1",
                              "BGN", "CA12", "LPIN1", "UBE2H"))
dev.off()


##----------------------------------------------------------##
#######----- Compare women 5 days with men 5 days -----#######
##----------------------------------------------------------##

### upload protein counts before NA imputation
data <- read.csv("./Area_data_42_samples_full_spec_lib.csv", sep = ",")

#metadata
experimental_design <- read_xlsx("./Sample_data_days.xlsx", sheet = 4)

#subset samples of interest (women 5 days and men 5 days)
data <- data[,c("Protein.Ids", experimental_design$label)]

data[data == 0] = NA

#filtering
data$na.count <- apply(data, 1, function(x) sum(is.na(x)))
data$na.count
nrow(data)

data <- subset(data, data$na.count < 5)
nrow(data)


#convert uniprot ids into gene names
proteins <- data$Protein.Ids 
symbols <- mapIds(org.Hs.eg.db, keys = proteins,
                  column = c('SYMBOL'),
                  keytype = 'UNIPROT')
proteins_symbol <- data.frame(proteins, symbols)

counter=1
for(protein in proteins_symbol$symbols){
  if (is.na(protein)) {
    proteins_symbol[counter,2] <- proteins_symbol[counter,1]
  }
  counter = counter+1
}

# change uniprot ids to protein symbols where possible
data$Protein <- proteins_symbol$symbols

#remove contaminants
data <- data[!grepl("contam_sp", data$Protein.Ids),]

#checking for duplicated proteins
data$Protein.Ids %>% duplicated() %>% any()

#prepare data for DEP analysis
data_unique <- make_unique(data, "Protein.Ids", "Protein", delim = ";")

#checking for duplicated proteins
data_unique$name %>% duplicated() %>% any()

#generate a SummarizedExperiment object using an experimental design
#and perform log2-transformation
ALM_columns <- grep("ALM", colnames(data_unique))
data_se <- make_se(data_unique, ALM_columns, experimental_design)

#plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)


#checking for missing values
data_filt <- filter_missval(data_se, thr = 1)
plot_frequency(data_filt)

#protein number
plot_numbers(data_se)
dev.off()

#normalize the data (if necessary)
data_norm <- normalize_vsn(data_filt)

#visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm) 

#plot a heatmap of proteins with missing values
plot_missval(data_filt)

#impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "knn")

#plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)


## Differential enrichment analysis  based on linear models and empherical Bayes statistics
#test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "woman", design_formula = formula(~ 0 + condition))  # specify the MSCs type for the "control" option

#denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05)

#save dep results as a file
dep_results <- get_results(dep)
colnames(dep_results)
dep_results %>% filter(significant) %>% nrow()  ## no significant DEPs

#write.csv(x = dep_results, file = "M5d_vs_W5d_prot_dep.csv") 

dep$replicate <- ""

## PCA
# tiff(file="./M5d_VS_W5d_PCA.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4, label = T) + scale_color_brewer(palette = "Set2")
dev.off()


#PLS-DA 
ordination.optimum.splsda <- splsda(t(data_norm@assays@data@listData[[1]]), experimental_design$condition, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

# tiff('PLSDA_M5d_vs_W5d_comp.tiff',
#      units="in",
#      width=10,
#      height=8,
#      res=300,
#      compression = 'lzw')
layout(matrix(c(1, 2)))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
dev.off()

# tiff('PLSDA_M5d_vs_W5d.tiff',
#      units="in",
#      width=10,
#      height=8,
#      res=300,
#      compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 2, size.axis = 1.0, size.xlabel = 1.2, size.ylabel = 1.2, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()

rownames(dep) <- dep@elementMetadata@listData[["Protein"]]



##----------------------------------------------------------##
######----- Compare women 10 days with men 10 days -----######
##----------------------------------------------------------##

### upload protein counts before NA imputation
data <- read.csv("./Area_data_42_samples_full_spec_lib.csv", sep = ",")

#metadata
experimental_design <- read_xlsx("./Sample_data_days.xlsx", sheet = 4)

#subset samples of interest (women 5 days and 10 days)
data <- data[,c("Protein.Ids", experimental_design$label)]

data[data == 0] = NA

#filtering
data$na.count <- apply(data, 1, function(x) sum(is.na(x)))
data$na.count
nrow(data)

data <- subset(data, data$na.count < 5)
nrow(data)


#convert uniprot ids into gene names
proteins <- data$Protein.Ids 
symbols <- mapIds(org.Hs.eg.db, keys = proteins,
                  column = c('SYMBOL'),
                  keytype = 'UNIPROT')
proteins_symbol <- data.frame(proteins, symbols)

counter=1
for(protein in proteins_symbol$symbols){
  if (is.na(protein)) {
    proteins_symbol[counter,2] <- proteins_symbol[counter,1]
  }
  counter = counter+1
}

# change uniprot ids to protein symbols where possible
data$Protein <- proteins_symbol$symbols

#remove contaminants
data <- data[!grepl("contam_sp", data$Protein.Ids),]

#checking for duplicated proteins
data$Protein.Ids %>% duplicated() %>% any()

#prepare data for DEP analysis
data_unique <- make_unique(data, "Protein.Ids", "Protein", delim = ";")

#checking for duplicated proteins
data_unique$name %>% duplicated() %>% any()

#generate a SummarizedExperiment object using an experimental design
#and perform log2-transformation
ALM_columns <- grep("ALM", colnames(data_unique))
data_se <- make_se(data_unique, ALM_columns, experimental_design)

#plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)


#checking for missing values
data_filt <- filter_missval(data_se, thr = 1)
plot_frequency(data_filt)

#protein number
plot_numbers(data_se)
dev.off()

#normalize the data (if necessary)
data_norm <- normalize_vsn(data_filt)

#visualize normalization by boxplots for all samples before and after normalization
plot_normalization(data_filt, data_norm) 

#plot a heatmap of proteins with missing values
plot_missval(data_filt)

#impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "knn")

#plot intensity distributions before and after imputation
plot_imputation(data_norm, data_imp)


## Differential enrichment analysis  based on linear models and empherical Bayes statistics
#test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "woman", design_formula = formula(~ 0 + condition))  # specify the MSCs type for the "control" option

#denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05)

#save dep results as a file
dep_results <- get_results(dep)
colnames(dep_results)
dep_results %>% filter(significant) %>% nrow()

#write.csv(x = dep_results, file = "M10d_vs_W10d_prot_dep.csv") 

dep$replicate <- ""

## PCA
# tiff(file="./M10d_vs_W10d_PCA.tiff",  # change file name
#     units = "in",
#     width = 5,
#     height = 3,
#     res=300)
plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4, label = T) + scale_color_brewer(palette = "Set2")
dev.off()


#PLS-DA 
ordination.optimum.splsda <- splsda(t(data_norm@assays@data@listData[[1]]), experimental_design$condition, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

# tiff('PLSDA_M10d_vs_W10d_comp.tiff',
#      units="in",
#      width=10,
#      height=8,
#      res=300,
#      compression = 'lzw')
layout(matrix(c(1, 2)))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
dev.off()

# tiff('PLSDA_M10d_vs_W10d.tiff',
#      units="in",
#      width=10,
#      height=8,
#      res=300,
#      compression = 'lzw')
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 2, size.axis = 1.0, size.xlabel = 1.2, size.ylabel = 1.2, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
dev.off()



####----------M & W Overlap ----------

men <- read.csv("./M10d_vs_M5d_prot_dep.csv", row.names = 1)
men_significant <- subset(men, men$X10d_vs_X5d_significant == TRUE)
men_proteins <- men_significant$ID

women <- read.csv("./W10d_vs_W5d_prot_dep.csv", row.names = 1)
women_significant <- subset(women, women$X10d_vs_X5d_significant == TRUE)
women_proteins <- women_significant$ID

intersect(men_proteins, women_proteins)

setdiff(men_proteins, women_proteins)
setdiff(women_proteins, men_proteins)
