# Initial steps, importing libraries

library(future) #simple and uniform way of evaluating R expressions
library(ggplot2) #plotting library 
library(bigmemory)#manipulation of big objects
library(Seurat) #seurat library 
sessionInfo() #collect information about the current R session 

#plan("multiprocess", workers=8) #multicore evaluation, multisession evaluation 
#workers=8 number of parallel futures that can be active at the same time before blocking 
#plan("multisession", workers=8)

#output directory to save the results
outdir <- '~/capsule/results/preprocessing_variable_features/'
dir.create(outdir)

## Load RDS files 
OVCA3.S <- readRDS("~/capsule/data/Rstudio_analysis/preprocessing_variable_features/OVCA3_S_UnprocessedCohort_SeuratObj_20191023.rds")
OVCA3.R <- readRDS("~/capsule/data/Rstudio_analysis/preprocessing_variable_features/OVCA3_R_UnprocessedCohort_SeuratObj_20191028.rds")

#Update Seurat objects for visualization 
OVCA3.S <- UpdateSeuratObject(object = OVCA3.S)
OVCA3.R <- UpdateSeuratObject(object = OVCA3.R)

dim(OVCA3.S)
#[1] 33538  1766
dim(OVCA3.R)
#[1] 33538   582

## Pair S and R data but also keep them separate
#OVCA3
OVCA3 <- merge(x=OVCA3.S, y=OVCA3.R, add.cell.id1 = "OVCA3_S", add.cell.id2 = "OVCA3_R", merge.data=TRUE, project = "OVCA3_Pair")
dim(OVCA3)
# 33538  2348
saveRDS(OVCA3, file = "~/capsule/results/preprocessing_variable_features/OVCA3_S+R_UnprocessedCohort_SeuratObj_2023.rds")

# STEP 2: Quality control 
##OVCA3 merge
#pattern chosen is anything that contains ^MT-
mito.genes <- grep(pattern = "^MT-", OVCA3@assays$RNA@counts@Dimnames[[1]], value = TRUE) 

#Calculate the percentage of mito genes. Get the sum of the columns of a matrix. 
percent.mito <- Matrix::colSums(x=GetAssayData(object=OVCA3, slot='counts')[mito.genes, ])/Matrix::colSums(x=GetAssayData(object=OVCA3, slot='counts'))

#Add a column to meta.data 
OVCA3[['percent.mito']] <- percent.mito

#Create a plot and save it as a pdf for the nFeature_RNA 
pdf(paste(outdir,"OVCA3_S+P_ViolinPlot_nFeatures_2023.pdf"))
vPlot <- VlnPlot(object = OVCA3, features = c("nFeature_RNA"))
print(vPlot)
dev.off() 

#Create a plot and save it as a pdf for the percent.mito activity 
pdf(paste(outdir,"OVCA3_S+P_ViolinPlot_perMito_2023.pdf"))
vPlot <- VlnPlot(object = OVCA3, features = c("percent.mito"))
print(vPlot)
dev.off() 
### Filter cells that have >5% mitochondrial counts 

#View all the plots at the same time 
pdf(paste(outdir,"OVCA3_S+P_ViolinPlot_ALL_2023.pdf"))
vPlotALL <- VlnPlot(OVCA3, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol=3)
print(vPlotALL)
dev.off()

#Create scatter plots to show a trend between two features 
#Create a plot for the RNA counts and the mito percentage
pdf(paste(outdir,"OVCA3_S+P_ScatterPlot_RNAcount-vs-mito_2023.pdf"))
fsPlot <- FeatureScatter(object=OVCA3, feature1 = "nCount_RNA", feature2 = "percent.mito")
print(fsPlot)
dev.off()

#Create a plot for the number of genes (count) vs the amount of RNA 
pdf(paste(outdir,"OVCA3_S+P_ScatterPlot_nGenes-vs-nRNA_2023.pdf"))
fsPlot <- FeatureScatter(object=OVCA3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(fsPlot)
dev.off()
## Curve up more Count more feature RNA, as counts get larger, features also get larger. The higher #of different genes the more cells   


## Descriptive stats:
summary(OVCA3@meta.data$nFeature_RNA) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#     45    3902    4488    4205    5009    8228

#Set up an upper and lower limit, na.rm=TRUE because the features column has some missing values 
## Upper limit is the mean of the feature values + 2 sd 
nFeatUpper <- mean(OVCA3@meta.data$nFeature_RNA, na.rm=TRUE) + 2*sd(OVCA3@meta.data$nFeature_RNA, na.rm=TRUE) ##  6208.011
nFeatLower <- 200 

summary(OVCA3@meta.data$percent.mito)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.00000 0.03749 0.04899 0.06257 0.06131 0.93613 
### assumption: cancer cells may have higher metabolic rate
##Set an upper limit with the mean and 2 sd, not use the 5%, we need to cut off a higher percentage anyways because cancer cells have a high metabolic rate 
perMitoUpper <- mean(OVCA3@meta.data$percent.mito, na.rm=TRUE) + 2*sd(OVCA3@meta.data$percent.mito, na.rm=TRUE) ## [1] 0.1471598

#Create more Vlnplots now adding lines for the limits established 
pdf(paste(outdir,"OVCA3_S+P_ViolinPlot_perMito_wCutoffs_2023.pdf"))
vPlot <- VlnPlot(object=OVCA3, features = c("percent.mito"))+geom_hline(yintercept=perMitoUpper, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()

pdf(paste(outdir,"OVCA3_S+P_ViolinPlot_nFeatures_wCutoffs_2023.pdf"))
vPlot <- VlnPlot(object=OVCA3, features = c("nFeature_RNA"))+geom_hline(yintercept=nFeatUpper, linetype="dashed", color="red", size=1)+geom_hline(yintercept=nFeatLower, linetype="dashed", color="red", size=1)
print(vPlot)
dev.off()

#Divide the data based on the limits established
## subset the paired dataset
OVCA3 <- subset(x=OVCA3, subset=nFeature_RNA > nFeatLower & nFeature_RNA < nFeatUpper & percent.mito < perMitoUpper)

dim(OVCA3)
# [1] 33538  2262
saveRDS(OVCA3, file = "~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_SeuratObj_2023.rds")

# STEP 3: Normalize
OVCA3 <- NormalizeData(object=OVCA3, normalization.method="LogNormalize", scale.factor=1e4)
dim(OVCA3)
# 33538  2262
saveRDS(OVCA3, file = "~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_SeuratObj_thruStep3_2023.rds")


# STEP 4: Detection of variable features across the single cells 
OVCA3 <- FindVariableFeatures(object=OVCA3, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf), do.plot=TRUE)
length(x = VariableFeatures(object=OVCA3))  
## number of variable features detected after filtering, normalization: 4287
dim(OVCA3) #3538  2262
saveRDS(OVCA3, "~/capsule/results/preprocessing_variable_features/OVCA3_QC-Subset-Cohort_SeuratObj_thruStep4_2023.rds")

# STEP 5: Scaling data and removing unwanted sources of variation 

OVCA3 <- ScaleData(object=OVCA3, vars.to.regress = c("nCount_RNA","percent.mito")) 
saveRDS(OVCA3, "~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_SeuratObj_thruStep5_2023.rds")

# STEP 6: Perform linear dimension reduction/ Determine PCA 
## OVCA3 
OVCA3 <- RunPCA(object=OVCA3, features=VariableFeatures(object=OVCA3), verbose = TRUE, npcs=200)

pdf(paste(outdir,"OVCA3_S+P_QC-Subset-Cohort_PCA-Dim-Loadings_2023.pdf"))
vdlPlot <- VizDimLoadings(object=OVCA3, dims = 1:2) #Only 2 principal components 
print(vdlPlot)
dev.off()
#Shows only two principal components and the genes that are being upregulated 

pdf(paste(outdir,"OVCA3_S+P_QC-Subset-Cohort_PCA-Dim-Plot_2023.pdf"))
dPlot <- DimPlot(object=OVCA3)
print(dPlot)
dev.off()
#Shows a PC plot, similar to UMAP, showing the correlation between 2 principal components, look at the genes expressed in the middle of both clusters 

OVCA3 <- ProjectDim(object=OVCA3)

pdf(paste(outdir,"OVCA3_S+P_QC-Subset-Cohort_PCA-Dim-Heatmaps-1thr12_2023.pdf"))
dhPlot <- DimHeatmap(object=OVCA3, dims = 1:9, cells = 500, balanced = TRUE)
print(dhPlot)
dev.off()
#Shows heatmaps for all the principal components 

saveRDS(OVCA3, "~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_SeuratObj_thruStep6_20238.rds")

# STEP 7: Find Clusters
library(data.table)
library(dplyr)
#library(rdetools) #might not be available 
library(clustree)

##OVCA3
OVCA3 <- FindNeighbors(object=OVCA3, dims = 1:50)
OVCA3 <- FindClusters(object=OVCA3, resolution = 0.6)
# Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck

#Number of nodes: 2262
#Number of edges: 97506

#Running Louvain algorithm...
#Maximum modularity in 10 random starts: 0.7501
#Number of communities: 5
#Elapsed time: 0 seconds

saveRDS(OVCA3, "~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_SeuratObj_thruStep7_Res0p6only_2023.rds")

file1 <- "~/capsule/results/preprocessing_variable_features/OVCA3OutputDataNames_res0p6_2023.csv";
file2 <- "~/capsule/results/preprocessing_variable_features/OVCA3OutputData_res0p6_2023.csv";

cat(OVCA3@meta.data$RNA_snn_res.0.6, file=file2, sep=",\n")

listNames <- OVCA3@meta.data$orig.ident
write.table(data.frame(listNames),
            row.names=FALSE,
            col.names = FALSE, 
            file = file1,
            sep=",")

mydat1 <- read.csv(file2)
mydat2 <- read.csv(file1)
fulldat <- cbind(mydat2[1],mydat1[1])
fulltab <- as.data.table(fulldat)
# name table columns
names(fulltab)[1] <- paste("UMI")
names(fulltab)[2] <- paste("cluster")
# group data based on clusters
group_by(fulltab, cluster)
# create a table counting unique UMIs/cells per cluster
tabPerClus <- fulltab %>% group_by(cluster) %>% count()
type <- sub("\\_.*","",fulltab$UMI)
fulltab <- cbind(fulltab, type)

divout <-  "~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_CellBreakdown_PerClusterPerType_Res0p6only_2023.csv"
A <- fulltab %>% group_by(cluster) %>% count(UMI)
write.csv(A, file=divout)

rm(file1,file2,mydat1,mydat2,fulldat,fulltab,tabPerClus,type,A)

# STEP 8: UMAP
##OVCA3 merge
OVCA3 <- RunUMAP(object=OVCA3, reduction = "pca", dims = 1:50)
pdf(paste(outdir,"OVCA3_S+P_QC-Subset-Cohort_UMAP_byCluster_Res0p6only_2023.pdf", width=12, height=8))
dPlot <- DimPlot(OVCA3, reduction="umap") + xlim(-15, 15) + ylim(-15, 15)
print(dPlot)
dev.off()
pdf(paste(outdir,"OVCA3_S+P_QC-Subset-Cohort_UMAP_clustersByType_Res0p6only_2023.pdf", width=12, height=8))
dPlot <- DimPlot(OVCA3, reduction="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
print(dPlot)
dev.off()
pdf(paste(outdir,"OVCA3_S+P_QC-Subset-Cohort_UMAP_clustersByType_noLegend_Res0p6only_2023.pdf", width=12, height=8))
dPlot <- DimPlot(OVCA3, reduction="umap", group.by='orig.ident') + theme(legend.position = "none") + xlim(-15, 15) + ylim(-15, 15)
print(dPlot)
dev.off()
postscript(file="~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_UMAP_bySampleID_Res0p6only_2023_forAdobe.eps", width=12, height=8, horizontal = TRUE, onefile = FALSE)
DimPlot(OVCA3, reduction="umap", group.by='orig.ident') + xlim(-15, 15) + ylim(-15, 15)
dev.off()
saveRDS(OVCA3, file = "~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_SeuratObj_thruStep8_Res0p6only_2023.rds")

#save embeddings
embeddings <- OVCA3@reductions[["umap"]]@cell.embeddings
clusters <- OVCA3@meta.data[["seurat_clusters"]]
emb_clus <- cbind(embeddings,clusters)
write.csv(emb_clus, sep=",","~/capsule/results/preprocessing_variable_features/OVCA3_UMAP_embeddings+clusters.csv")

#STEP DEG: Differential Expressed Genes
## OVCA3
OVCA3.markers <- FindAllMarkers(OVCA3, min.pct = 0.25, logfc.threshold = 0.25)

saveRDS(OVCA3.markers, file = "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_ALL_2023.rds" )
dim(OVCA3.markers)
# [1] 9126  7
# Save
write.csv(OVCA3.markers, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_ALL_2023.csv")

#split object to check sensitive and persistent genes 
batches <- SplitObject(OVCA3, split.by = "orig.ident")
OVCA3.S <- batches[[1]]
OVCA3.R <- batches[[2]]
##find markers for sensitive and persistent 
OVCA3.S.markers <- FindAllMarkers(OVCA3.S, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(OVCA3.S.markers, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3.S_ALL_2023.csv")

OVCA3.R.markers <- FindAllMarkers(OVCA3.R, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(OVCA3.R.markers, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3.R_ALL_2023.csv")

## clusters: 0 - 4 
OVCA3clus0markers <- OVCA3.markers %>% group_by(cluster) %>% filter(cluster == 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus0markers, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-0_2023.csv")
###Only UP genes 
OVCA3clus0markers.UP <-OVCA3clus0markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC > 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus0markers.UP, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-0_UP_2023.csv")
###Only DOWN genes 
OVCA3clus0markers.DOWN <-OVCA3clus0markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC < 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus0markers.DOWN, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-0_DOWN_2023.csv")


OVCA3clus1markers <- OVCA3.markers %>% group_by(cluster) %>% filter(cluster == 1) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus1markers, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-1_2023.csv")
###Only UP genes 
OVCA3clus1markers.UP <-OVCA3clus1markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC > 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus1markers.UP, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-1_UP_2023.csv")
###Only DOWN genes 
OVCA3clus1markers.DOWN <-OVCA3clus1markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC < 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus1markers.DOWN, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-1_DOWN_2023.csv")


OVCA3clus2markers <- OVCA3.markers %>% group_by(cluster) %>% filter(cluster == 2) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus2markers, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-2_2023.csv")
###Only UP genes 
OVCA3clus2markers.UP <-OVCA3clus2markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC > 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus2markers.UP, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-2_UP_2023.csv")
###Only DOWN genes 
OVCA3clus2markers.DOWN <-OVCA3clus2markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC < 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus2markers.DOWN, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-2_DOWN_2023.csv")


OVCA3clus3markers <- OVCA3.markers %>% group_by(cluster) %>% filter(cluster == 3) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus3markers, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-3_2023.csv")
###Only UP genes 
OVCA3clus3markers.UP <-OVCA3clus3markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC > 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus3markers.UP, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-3_UP_2023.csv")
###Only DOWN genes 
OVCA3clus3markers.DOWN <-OVCA3clus3markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC < 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus3markers.DOWN, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-3_DOWN_2023.csv")


OVCA3clus4markers <- OVCA3.markers %>% group_by(cluster) %>% filter(cluster == 4) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus4markers, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-4_2023.csv")
###Only UP genes 
OVCA3clus4markers.UP <-OVCA3clus4markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC > 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus4markers.UP, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-4_UP_2023.csv")
###Only DOWN genes 
OVCA3clus4markers.DOWN <-OVCA3clus4markers %>% group_by(avg_log2FC) %>% filter(avg_log2FC < 0) %>% select(avg_log2FC,gene)
write.csv(OVCA3clus4markers.DOWN, "~/capsule/results/preprocessing_variable_features/DiffExpGeneList_OVCA3_Cluster-4_DOWN_2023.csv")

