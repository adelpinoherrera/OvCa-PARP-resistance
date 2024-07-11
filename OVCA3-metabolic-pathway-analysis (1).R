
library(future)
library(ggplot2)
library(bigmemory)
library(Seurat)
library(data.table)
library(dplyr)
#library(rdetools)
library(clustree)
sessionInfo()

#plan("multisession", workers=16)

#options(future.globals.maxSize=134217728000) ### for 128GB * 1024*1024

date = "_2023"
#output directory
outdir <- '~/capsule/results/metabolic_pathway_analysis/'
dir.create(outdir)

#load OVCA3 S+P object 
OVCA3 <- readRDS("~/capsule/results/preprocessing_variable_features/OVCA3_S+R_QC-Subset-Cohort_SeuratObj_thruStep8_Res0p6only_2023.rds")
#OVCA3 <- readRDS("~capsule/data/Rstudio_analysis/metabolic_pathway_analysis/OVCA3_S+R_QC-Subset-Cohort_SeuratObj_thruStep8_Res0p6only_2023.rds")
#OVCA3
names(OVCA3[[]])

#"orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mito"    "RNA_snn_res.0.6" "seurat_clusters"
OVCA3$Cluster <- OVCA3$RNA_snn_res.0.6
OVCA3$name <- OVCA3$orig.ident
names(OVCA3[[]])
#"orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mito"    "RNA_snn_res.0.6" "seurat_clusters" "Cluster"         "name"  


### Labeling as sensitive and resistant, create structure of meta data based on name, to be filled with new annotations
nameListMtx <- as.matrix(OVCA3$name)
tmp1 <- as.matrix(OVCA3$name)
numCells <- 2262
for (val in 1:numCells) {
  ### check name
  tmpname <- OVCA3[[]][val,'name']
  ## from name set CMML or normal
  if (isTRUE(tmpname == "OVCA3_R")){
    tmp1[val] <- "resistant"
  } else {
    tmp1[val] <- "sensitive"
  }
}

OVCA3Grpd <- OVCA3
#OVCA3Grpd.1$samptype <- tmp1

names(OVCA3Grpd[[]])
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mito"    "RNA_snn_res.0.6" "seurat_clusters" "S.Score"         "G2M.Score"      
# [9] "Phase"           "old.ident"       "Cluster"         "name" 



OVCA3_S_cells <- WhichCells(OVCA3, expression = (orig.ident=="OVCA3_S"))
OVCA3_R_cells <- WhichCells(OVCA3, expression = (orig.ident=="OVCA3_R"))


celllinesonly <- c(OVCA3_S_cells, OVCA3_R_cells)

OVCA3CellLineGrpd <- subset(x=OVCA3Grpd, cells=celllinesonly)

OVCA3CellLineGrpd@active.ident <- OVCA3CellLineGrpd$Cluster

Sonly <- c(OVCA3_S_cells)
Ronly <- c(OVCA3_R_cells)

OVCA3pair <- subset(x=OVCA3CellLineGrpd, cells=c(OVCA3_S_cells,OVCA3_R_cells))

##Metabolic Pathway score analysis for 42 different pathways
#1
OXPHOSgenes <- c('ATP12A', 'ATP4A', 'ATP4B', 'ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D', 'ATP5E', 'ATP5F1', 'ATP5G1', 'ATP5G1P5', 'ATP5G2', 'ATP5G3', 'ATP5H', 'ATP5I', 'ATP5J', 'ATP5J2', 'ATP5L', 'ATP5O', 'ATP6', 'ATP6AP1', 'ATP6V0A1', 'ATP6V0A2', 'ATP6V0A4', 'ATP6V0B', 'ATP6V0C', 'ATP6V0D1', 'ATP6V0D2', 'ATP6V0E1', 'ATP6V0E2', 'ATP6V1A', 'ATP6V1B1', 'ATP6V1B2', 'ATP6V1C1', 'ATP6V1C2', 'ATP6V1D', 'ATP6V1E1', 'ATP6V1E2', 'ATP6V1F', 'ATP6V1G1', 'ATP6V1G2', 'ATP6V1G3', 'ATP6V1H', 'ATP8', 'COX1', 'COX10', 'COX11', 'COX15', 'COX17', 'COX2', 'COX3', 'COX4I1', 'COX4I2', 'COX5A', 'COX5B', 'COX6A1', 'COX6A2', 'COX6B1', 'COX6B2', 'COX6C', 'COX6CP3', 'COX7A1', 'COX7A2', 'COX7A2L', 'COX7B', 'COX7B2', 'COX7C', 'COX8A', 'COX8C', 'CYC1', 'CYTB', 'LHPP', 'LOC100133737', 'LOC642502', 'LOC644310', 'LOC727947', 'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA4L2', 'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFB8', 'NDUFB9', 'NDUFC1', 'NDUFC2', 'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5', 'NDUFS6', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFV3', 'PPA1', 'PPA2', 'SDHA', 'SDHB', 'SDHC', 'SDHD', 'TCIRG1', 'UQCR10', 'UQCR11', 'UQCRB', 'UQCRC1', 'UQCRC2', 'UQCRFS1', 'UQCRH', 'UQCRHL', 'UQCRQ')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "OXPHOS", search=TRUE)

#2
glycolysisGenes <- c('ACSS1', 'ACSS2', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7', 'AKR1A1', 'ALDH1A3', 'ALDH1B1', 'ALDH2', 'ALDH3A1', 'ALDH3A2', 'ALDH3B1', 'ALDH3B2', 'ALDH7A1', 'ALDH9A1', 'ALDOA', 'ALDOB', 'ALDOC', 'BPGM', 'DLAT', 'DLD', 'ENO1', 'ENO2', 'ENO3', 'FBP1', 'FBP2', 'G6PC', 'G6PC2', 'GALM', 'GAPDH', 'GCK', 'GPI', 'HK1', 'HK2', 'HK3', 'LDHA', 'LDHAL6A', 'LDHAL6B', 'LDHB', 'LDHC', 'PCK1', 'PCK2', 'PDHA1', 'PDHA2', 'PDHB', 'PFKL', 'PFKM', 'PFKP', 'PGAM1', 'PGAM2', 'PGAM4', 'PGK1', 'PGK2', 'PGM1', 'PGM2', 'PKLR', 'PKM2', 'TPI1')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "glycolysis", search=TRUE)

#3
citricAcidGenes <- c('ACLY', 'ACO1', 'ACO2', 'CS', 'DLAT', 'DLD', 'DLST', 'FH', 'IDH1', 'IDH2', 'IDH3A', 'IDH3B', 'IDH3G', 'LOC283398', 'LOC642502', 'MDH1', 'MDH2', 'OGDH', 'OGDHL', 'PC', 'PCK1', 'PCK2', 'PDHA1', 'PDHA2', 'PDHB', 'SDHA', 'SDHB', 'SDHC', 'SDHD', 'SUCLA2', 'SUCLG1', 'SUCLG2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "citricAcid", search=TRUE)

#4
ascorbAldarMet <- c('ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH7A1', 'ALDH9A1', 'MIOX', 'UGDH', 'UGT1A1', 'UGT1A10', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT2A1', 'UGT2A3', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28', 'UGT2B4', 'UGT2B7')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "ascorbateAldarate", search=TRUE)

#5
chondroitinSulfateBios <- c('B3GALT6', 'B3GAT1', 'B3GAT2', 'B3GAT3', 'B4GALT7', 'CHPF', 'CHPF2', 'CHST11', 'CHST12', 'CHST13', 'CHST14', 'CHST15', 'CHST3', 'CHST7', 'CHSY1', 'CHSY3', 'CSGALNACT1', 'CSGALNACT2', 'DSE', 'UST', 'XYLT1', 'XYLT2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "chondroitinSulfate", search=TRUE)

#6
fructMannMet <- c('AKR1B1', 'AKR1B10', 'ALDOA', 'ALDOB', 'ALDOC', 'FBP1', 'FBP2', 'FPGT', 'FUK', 'GMDS', 'GMPPA', 'GMPPB', 'HK1', 'HK2', 'HK3', 'KHK', 'MPI', 'MTMR1', 'MTMR2', 'MTMR6', 'MTMR7', 'PFKFB1', 'PFKFB2', 'PFKFB3', 'PFKFB4', 'PFKL', 'PFKM', 'PFKP', 'PHPT1', 'PMM1', 'PMM2', 'SORD', 'TPI1', 'TSTA3')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "fructoseMannose", search=TRUE)

#7
galactose <- c('AKR1B1', 'B4GALT1', 'B4GALT2', 'G6PC', 'G6PC2', 'GAA', 'GALE', 'GALK1', 'GALK2', 'GALT', 'GANC', 'GCK', 'GLA', 'GLB1', 'HK1', 'HK2', 'HK3', 'LALBA', 'LCT', 'MGAM', 'PFKL', 'PFKM', 'PFKP', 'PGM1', 'PGM2', 'UGP2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "galactose", search=TRUE)

#8
keratanSulfate <- c('B3GNT1', 'B3GNT2', 'B3GNT7', 'B4GALT1', 'B4GALT2', 'B4GALT3', 'B4GALT4', 'CHST1', 'CHST2', 'CHST4', 'CHST6', 'FUT8', 'ST3GAL1', 'ST3GAL2', 'ST3GAL3')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "keratanSulfate", search=TRUE)

#9
heparanSulfate <- c('B3GALT6', 'B3GAT1', 'B3GAT2', 'B3GAT3', 'B4GALT7', 'EXT1', 'EXT2', 'EXTL1', 'EXTL2', 'EXTL3', 'GLCE', 'HS2ST1', 'HS3ST1', 'HS3ST2', 'HS3ST3A1', 'HS3ST3B1', 'HS3ST5', 'HS6ST1', 'HS6ST2', 'HS6ST3', 'NDST1', 'NDST2', 'NDST3', 'NDST4', 'XYLT1', 'XYLT2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "heparanSulfate", search=TRUE)

#10
pentosePhosphate <- c('ALDOA', 'ALDOB', 'ALDOC', 'DERA', 'FBP1', 'FBP2', 'G6PD', 'GPI', 'H6PD', 'LOC729020', 'PFKL', 'PFKM', 'PFKP', 'PGD', 'PGLS', 'PGM1', 'PGM2', 'PRPS1', 'PRPS1L1', 'PRPS2', 'RBKS', 'RPE', 'RPIA', 'TALDO1', 'TKT', 'TKTL1', 'TKTL2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "pentosePhosphate", search=TRUE)

#11
pyruvateMet <- c('ACACA', 'ACACB', 'ACAT1', 'ACAT2', 'ACOT12', 'ACSS1', 'ACSS2', 'ACYP1', 'ACYP2', 'AKR1B1', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH7A1', 'ALDH9A1', 'DLAT', 'DLD', 'GLO1', 'GRHPR', 'HAGH', 'HAGHL', 'LDHA', 'LDHAL6A', 'LDHAL6B', 'LDHB', 'LDHC', 'LDHD', 'MDH1', 'MDH2', 'ME1', 'ME2', 'ME3', 'PC', 'PCK1', 'PCK2', 'PDHA1', 'PDHA2', 'PDHB', 'PKLR', 'PKM2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "pyruvateMet", search=TRUE)

#12
starchSucrose <- c('AGL', 'AMY1A', 'AMY1B', 'AMY1C', 'AMY2A', 'AMY2B', 'ENPP1', 'ENPP3', 'G6PC', 'G6PC2', 'GAA', 'GANC', 'GBA3', 'GBE1', 'GCK', 'GPI', 'GUSB', 'GYS1', 'GYS2', 'HK1', 'HK2', 'HK3', 'MGAM', 'PGM1', 'PGM2', 'PGM2L1', 'PYGB', 'PYGL', 'PYGM', 'SI', 'TREH', 'UGDH', 'UGP2', 'UGT1A1', 'UGT1A10', 'UGT1A3', 'UGT1A4', 'UGT1A5', 'UGT1A6', 'UGT1A7', 'UGT1A8', 'UGT1A9', 'UGT2A1', 'UGT2A3', 'UGT2B10', 'UGT2B11', 'UGT2B15', 'UGT2B17', 'UGT2B28', 'UGT2B4', 'UGT2B7', 'UXS1')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "starchSucrose", search=TRUE)

#13
etherLipidMet <- c('AGPS', 'CHPT1', 'ENPP2', 'ENPP6', 'JMJD7-PLA2G4B', 'LPCAT1', 'LPCAT2', 'LPCAT4', 'PAFAH1B1', 'PAFAH1B2', 'PAFAH1B3', 'PAFAH2', 'PLA2G10', 'PLA2G12A', 'PLA2G12B', 'PLA2G1B', 'PLA2G2A', 'PLA2G2C', 'PLA2G2D', 'PLA2G2E', 'PLA2G2F', 'PLA2G3', 'PLA2G4A', 'PLA2G4B', 'PLA2G4E', 'PLA2G5', 'PLA2G6', 'PLA2G7', 'PLD1', 'PLD2', 'PPAP2A', 'PPAP2B', 'PPAP2C')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "etherLipidMet", search=TRUE)

#14
arachidonicAcidMet <- c('AKR1C3', 'ALOX12', 'ALOX12B', 'ALOX15', 'ALOX15B', 'ALOX5', 'CBR1', 'CBR3', 'CYP2B6', 'CYP2C18', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP2E1', 'CYP2J2', 'CYP2U1', 'CYP4A11', 'CYP4A22', 'CYP4F2', 'CYP4F3', 'EPHX2', 'GGT1', 'GGT5', 'GGT6', 'GGT7', 'GPX1', 'GPX2', 'GPX3', 'GPX4', 'GPX5', 'GPX6', 'GPX7', 'HPGDS', 'JMJD7-PLA2G4B', 'LTA4H', 'LTC4S', 'PLA2G10', 'PLA2G12A', 'PLA2G12B', 'PLA2G1B', 'PLA2G2A', 'PLA2G2C', 'PLA2G2D', 'PLA2G2E', 'PLA2G2F', 'PLA2G3', 'PLA2G4A', 'PLA2G4B', 'PLA2G4E', 'PLA2G5', 'PLA2G6', 'PTGDS', 'PTGES', 'PTGES2', 'PTGIS', 'PTGS1', 'PTGS2', 'TBXAS1')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "arachidonicAcidMet", search=TRUE)

#15
fattyAcidMet <- c('ACAA1', 'ACAA2', 'ACADL', 'ACADM', 'ACADS', 'ACADSB', 'ACADVL', 'ACAT1', 'ACAT2', 'ACOX1', 'ACOX3', 'ACSL1', 'ACSL3', 'ACSL4', 'ACSL5', 'ACSL6', 'ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH7A1', 'ALDH9A1', 'CPT1A', 'CPT1B', 'CPT1C', 'CPT2', 'CYP4A11', 'CYP4A22', 'ECHS1', 'ECI1', 'ECI2', 'EHHADH', 'GCDH', 'HADH', 'HADHA', 'HADHB')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "fattyAcidMet", search=TRUE)

#16
glycerolipidMet <- c('AGK', 'AGPAT1', 'AGPAT2', 'AGPAT3', 'AGPAT4', 'AGPAT6', 'AGPAT9', 'AKR1A1', 'AKR1B1', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH7A1', 'ALDH9A1', 'AWAT2', 'CEL', 'DAK', 'DGAT1', 'DGAT2', 'DGKA', 'DGKB', 'DGKD', 'DGKE', 'DGKG', 'DGKH', 'DGKI', 'DGKQ', 'DGKZ', 'GK', 'GK2', 'GLA', 'GLYCTK', 'GPAM', 'GPAT2', 'LCLAT1', 'LIPC', 'LIPF', 'LIPG', 'LPL', 'MBOAT1', 'MBOAT2', 'MGLL', 'PNLIP', 'PNLIPRP1', 'PNLIPRP2', 'PNPLA3', 'PPAP2A', 'PPAP2B', 'PPAP2C')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "glycerolipidMet", search=TRUE)

#17
glycerophospholipidMet <- c('ACHE', 'AGPAT1', 'AGPAT2', 'AGPAT3', 'AGPAT4', 'AGPAT6', 'AGPAT9', 'C17orf48', 'CDIPT', 'CDS1', 'CDS2', 'CHAT', 'CHKA', 'CHKB', 'CHPT1', 'CRLS1', 'DGKA', 'DGKB', 'DGKD', 'DGKE', 'DGKG', 'DGKH', 'DGKI', 'DGKQ', 'DGKZ', 'ETNK1', 'ETNK2', 'GNPAT', 'GPAM', 'GPAT2', 'GPD1', 'GPD1L', 'GPD2', 'JMJD7-PLA2G4B', 'LCAT', 'LCLAT1', 'LPCAT1', 'LPCAT2', 'LPCAT3', 'LPCAT4', 'LPGAT1', 'LYPLA1', 'LYPLA2', 'MBOAT1', 'MBOAT2', 'MBOAT7', 'PCYT1A', 'PCYT1B', 'PCYT2', 'PEMT', 'PGS1', 'PHOSPHO1', 'PISD', 'PLA2G10', 'PLA2G12A', 'PLA2G12B', 'PLA2G15', 'PLA2G1B', 'PLA2G2A', 'PLA2G2C', 'PLA2G2D', 'PLA2G2E', 'PLA2G2F', 'PLA2G3', 'PLA2G4A', 'PLA2G4B', 'PLA2G4E', 'PLA2G5', 'PLA2G6', 'PLD1', 'PLD2', 'PPAP2A', 'PPAP2B', 'PPAP2C', 'PTDSS1', 'PTDSS2', 'TAZ')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "glycerophospholipidMet", search=TRUE)

#18
glyoxylateDicarbMet <- c('ACO1', 'ACO2', 'AFMID', 'CS', 'GLYCTK', 'GRHPR', 'HAO1', 'HAO2', 'HYI', 'MDH1', 'MDH2', 'MTHFD1', 'MTHFD1L', 'MTHFD2', 'MTHFD2L', 'PGP')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "glyoxylateDicarbMet", search=TRUE)

#19
linoleicAcidMet <- c('AKR1B10', 'ALOX15', 'CYP1A2', 'CYP2C18', 'CYP2C19', 'CYP2C8', 'CYP2C9', 'CYP2E1', 'CYP2J2', 'CYP3A4', 'CYP3A43', 'CYP3A5', 'CYP3A7', 'JMJD7-PLA2G4B', 'PLA2G10', 'PLA2G12A', 'PLA2G12B', 'PLA2G1B', 'PLA2G2A', 'PLA2G2C', 'PLA2G2D', 'PLA2G2E', 'PLA2G2F', 'PLA2G3', 'PLA2G4A', 'PLA2G4B', 'PLA2G4E', 'PLA2G5', 'PLA2G6')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "linoleicAcidMet", search=TRUE)

#20
sphingolipidMet <- c('ACER1', 'ACER2', 'ACER3', 'ARSA', 'ASAH1', 'ASAH2', 'ASAH2C', 'B4GALT6', 'CERK', 'DEGS1', 'DEGS2', 'ENPP7', 'GAL3ST1', 'GALC', 'GBA', 'GLA', 'GLB1', 'KDSR', 'NEU1', 'NEU2', 'NEU3', 'NEU4', 'PPAP2A', 'PPAP2B', 'PPAP2C', 'SGMS1', 'SGMS2', 'SGPL1', 'SGPP1', 'SGPP2', 'SMPD1', 'SMPD2', 'SMPD3', 'SMPD4', 'SPHK1', 'SPHK2', 'SPTLC1', 'SPTLC2', 'UGCG', 'UGT8')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "sphingolipidMet", search=TRUE)

#21
propanoateMet <- c('ABAT', 'ACACA', 'ACACB', 'ACADM', 'ACAT1', 'ACAT2', 'ACSS1', 'ACSS2', 'ACSS3', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH6A1', 'ALDH7A1', 'ALDH9A1', 'ECHS1', 'EHHADH', 'HADHA', 'HIBCH', 'LDHA', 'LDHAL6A', 'LDHAL6B', 'LDHB', 'LDHC', 'LOC283398', 'MCEE', 'MLYCD', 'MUT', 'PCCA', 'PCCB', 'SUCLA2', 'SUCLG1', 'SUCLG2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "propanoateMet", search=TRUE)

#22
pyrimidineMet <- c('AK3', 'CAD', 'CANT1', 'CDA', 'CMPK1', 'CMPK2', 'CTPS', 'CTPS2', 'DCK', 'DCTD', 'DHODH', 'DPYD', 'DPYS', 'DTYMK', 'DUT', 'ENTPD1', 'ENTPD3', 'ENTPD4', 'ENTPD5', 'ENTPD6', 'ENTPD8', 'ITPA', 'NME1', 'NME1-NME2', 'NME2', 'NME3', 'NME4', 'NME5', 'NME6', 'NME7', 'NT5C', 'NT5C1A', 'NT5C1B', 'NT5C2', 'NT5C3', 'NT5E', 'NT5M', 'NUDT2', 'PNP', 'PNPT1', 'POLA1', 'POLA2', 'POLD1', 'POLD2', 'POLD3', 'POLD4', 'POLE', 'POLE2', 'POLE3', 'POLE4', 'POLR1A', 'POLR1B', 'POLR1C', 'POLR1D', 'POLR1E', 'POLR2A', 'POLR2B', 'POLR2C', 'POLR2D', 'POLR2E', 'POLR2F', 'POLR2G', 'POLR2H', 'POLR2I', 'POLR2J', 'POLR2J2', 'POLR2J3', 'POLR2K', 'POLR2L', 'POLR3A', 'POLR3B', 'POLR3C', 'POLR3D', 'POLR3F', 'POLR3G', 'POLR3GL', 'POLR3H', 'POLR3K', 'PRIM1', 'PRIM2', 'RRM1', 'RRM2', 'RRM2B', 'TK1', 'TK2', 'TXNRD1', 'TXNRD2', 'TYMP', 'TYMS', 'UCK1', 'UCK2', 'UCKL1', 'UMPS', 'UPB1', 'UPP1', 'UPP2', 'UPRT', 'ZNRD1')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(OXPHOSgenes), name = "pyrimidineMet", search=TRUE)

#23
purineMet <- c('ADA', 'ADCY1', 'ADCY10', 'ADCY2', 'ADCY3', 'ADCY4', 'ADCY5', 'ADCY6', 'ADCY7', 'ADCY8', 'ADCY9', 'ADK', 'ADSL', 'ADSS', 'ADSSL1', 'AK1', 'AK2', 'AK4', 'AK5', 'AK7', 'ALLC', 'AMPD1', 'AMPD2', 'AMPD3', 'APRT', 'ATIC', 'C17orf48', 'CANT1', 'DCK', 'DGUOK', 'ENPP1', 'ENPP3', 'ENTPD1', 'ENTPD2', 'ENTPD3', 'ENTPD4', 'ENTPD5', 'ENTPD6', 'ENTPD8', 'FHIT', 'GART', 'GDA', 'GMPR', 'GMPR2', 'GMPS', 'GUCY1A2', 'GUCY1A3', 'GUCY1B3', 'GUCY2C', 'GUCY2D', 'GUCY2F', 'GUK1', 'HPRT1', 'IMPDH1', 'IMPDH2', 'ITPA', 'NME1', 'NME1-NME2', 'NME2', 'NME3', 'NME4', 'NME5', 'NME6', 'NME7', 'NPR1', 'NPR2', 'NT5C', 'NT5C1A', 'NT5C1B', 'NT5C2', 'NT5C3', 'NT5E', 'NT5M', 'NUDT2', 'NUDT5', 'NUDT9', 'PAICS', 'PAPSS1', 'PAPSS2', 'PDE10A', 'PDE11A', 'PDE1A', 'PDE1B', 'PDE1C', 'PDE2A', 'PDE3A', 'PDE3B', 'PDE4A', 'PDE4B', 'PDE4C', 'PDE4D', 'PDE5A', 'PDE6A', 'PDE6B', 'PDE6C', 'PDE6D', 'PDE6G', 'PDE6H', 'PDE7A', 'PDE7B', 'PDE8A', 'PDE8B', 'PDE9A', 'PFAS', 'PKLR', 'PKM2', 'PNP', 'PNPT1', 'POLA1', 'POLA2', 'POLD1', 'POLD2', 'POLD3', 'POLD4', 'POLE', 'POLE2', 'POLE3', 'POLE4', 'POLR1A', 'POLR1B', 'POLR1C', 'POLR1D', 'POLR1E', 'POLR2A', 'POLR2B', 'POLR2C', 'POLR2D', 'POLR2E', 'POLR2F', 'POLR2G', 'POLR2H', 'POLR2I', 'POLR2J', 'POLR2J2', 'POLR2J3', 'POLR2K', 'POLR2L', 'POLR3A', 'POLR3B', 'POLR3C', 'POLR3D', 'POLR3F', 'POLR3G', 'POLR3GL', 'POLR3H', 'POLR3K', 'PPAT', 'PRHOXNB', 'PRIM1', 'PRIM2', 'PRPS1', 'PRPS1L1', 'PRPS2', 'PRUNE', 'RRM1', 'RRM2', 'RRM2B', 'XDH', 'ZNRD1')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(purineMet), name = "purineMet", search=TRUE)

#24
alaniAspratGlutamMet <- c('ABAT', 'ACY3', 'ADSL', 'ADSS', 'ADSSL1', 'AGXT', 'AGXT2', 'ALDH4A1', 'ALDH5A1', 'ASL', 'ASNS', 'ASPA', 'ASS1', 'CAD', 'CPS1', 'DDO', 'GAD1', 'GAD2', 'GFPT1', 'GFPT2', 'GLS', 'GLS2', 'GLUD1', 'GLUD2', 'GLUL', 'GOT1', 'GOT2', 'GPT', 'GPT2', 'IL4I1', 'NIT2', 'PPAT')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(alaniAspratGlutamMet), name = "alaniAspratGlutamMet", search=TRUE)

#25
arginineProlineMet <- c('ABP1', 'ACY1', 'ADC', 'AGMAT', 'ALDH18A1', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH4A1', 'ALDH7A1', 'ALDH9A1', 'AMD1', 'ARG1', 'ARG2', 'ASL', 'ASS1', 'CKB', 'CKM', 'CKMT1A', 'CKMT1B', 'CKMT2', 'CPS1', 'DAO', 'GAMT', 'GATM', 'GLS', 'GLS2', 'GLUD1', 'GLUD2', 'GLUL', 'GOT1', 'GOT2', 'LAP3', 'MAOA', 'MAOB', 'NAGS', 'NOS1', 'NOS2', 'NOS3', 'OAT', 'ODC1', 'OTC', 'P4HA1', 'P4HA2', 'P4HA3', 'PRODH', 'PRODH2', 'PYCR1', 'PYCR2', 'PYCRL', 'SAT1', 'SAT2', 'SMS', 'SRM')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(arginineProlineMet), name = "arginineProlineMet", search=TRUE)

#26
butanoateMet <- c('AACS', 'ABAT', 'ACADS', 'ACAT1', 'ACAT2', 'ACSM1', 'ACSM2A', 'ACSM3', 'ACSM4', 'ACSM5', 'AKR1B10', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH5A1', 'ALDH7A1', 'ALDH9A1', 'BDH1', 'BDH2', 'ECHS1', 'EHHADH', 'GAD1', 'GAD2', 'HADH', 'HADHA', 'HMGCL', 'HMGCS1', 'HMGCS2', 'L2HGDH', 'OXCT1', 'OXCT2', 'PDHA1', 'PDHA2', 'PDHB')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(butanoateMet), name = "butanoateMet", search=TRUE)

#27
glutathioneMet <- c('ANPEP', 'G6PD', 'GCLC', 'GCLM', 'GGCT', 'GGT1', 'GGT5', 'GGT6', 'GGT7', 'GPX1', 'GPX2', 'GPX3', 'GPX4', 'GPX5', 'GPX6', 'GPX7', 'GSR', 'GSS', 'GSTA1', 'GSTA2', 'GSTA3', 'GSTA4', 'GSTA5', 'GSTK1', 'GSTM1', 'GSTM2', 'GSTM3', 'GSTM4', 'GSTM5', 'GSTO1', 'GSTO2', 'GSTP1', 'GSTT1', 'GSTT2', 'GSTZ1', 'IDH1', 'IDH2', 'LAP3', 'MGST1', 'MGST2', 'MGST3', 'ODC1', 'OPLAH', 'PGD', 'RRM1', 'RRM2', 'RRM2B', 'SMS', 'SRM', 'TXNDC12')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(glutathioneMet), name = "glutathioneMet", search=TRUE)

#28
glycineSerineThreonMet <- c('AGXT', 'AGXT2', 'ALAS1', 'ALAS2', 'AMT', 'AOC2', 'AOC3', 'BHMT', 'CBS', 'CHDH', 'CTH', 'DAO', 'DLD', 'DMGDH', 'GAMT', 'GATM', 'GCAT', 'GLDC', 'GLYCTK', 'GNMT', 'MAOA', 'MAOB', 'PHGDH', 'PIPOX', 'PSAT1', 'PSPH', 'SARDH', 'SDS', 'SHMT1', 'SHMT2', 'SRR')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(glycineSerineThreonMet), name = "glycineSerineThreonMet", search=TRUE)

#29
histidineMet <- c('ABP1', 'ACY3', 'ALDH1A3', 'ALDH1B1', 'ALDH2', 'ALDH3A1', 'ALDH3A2', 'ALDH3B1', 'ALDH3B2', 'ALDH7A1', 'ALDH9A1', 'AMDHD1', 'ASPA', 'CNDP1', 'DDC', 'FTCD', 'HAL', 'HDC', 'HEMK1', 'HNMT', 'LCMT1', 'LCMT2', 'MAOA', 'MAOB', 'METTL2B', 'METTL6', 'TRMT11', 'UROC1', 'WBSCR22')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(histidineMet), name = "histidineMet", search=TRUE)

#30
lysineDeg <- c('AADAT', 'AASDH', 'AASDHPPT', 'AASS', 'ACAT1', 'ACAT2', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH7A1', 'ALDH9A1', 'ASH1L', 'BBOX1', 'DLST', 'DOT1L', 'ECHS1', 'EHHADH', 'EHMT1', 'EHMT2', 'GCDH', 'HADH', 'HADHA', 'NSD1', 'OGDH', 'OGDHL', 'PIPOX', 'PLOD1', 'PLOD2', 'PLOD3', 'SETD1A', 'SETD1B', 'SETD2', 'SETD7', 'SETD8', 'SETDB1', 'SETDB2', 'SETMAR', 'SUV39H1', 'SUV39H2', 'SUV420H1', 'SUV420H2', 'TMLHE', 'WHSC1', 'WHSC1L1')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(lysineDeg), name = "lysineDeg", search=TRUE)

#31
cysteineMethionineMet <- c('ADI1', 'AHCY', 'AHCYL1', 'AHCYL2', 'AMD1', 'APIP', 'BHMT', 'CBS', 'CDO1', 'CTH', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DNMT3L', 'ENOPH1', 'GOT1', 'GOT2', 'IL4I1', 'LDHA', 'LDHAL6A', 'LDHAL6B', 'LDHB', 'LDHC', 'MAT1A', 'MAT2A', 'MAT2B', 'MPST', 'MTAP', 'MTR', 'SDS', 'SMS', 'SRM', 'TAT', 'TRDMT1')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(cysteineMethionineMet), name = "cysteineMethionineMet", search=TRUE)

#32
phenylalanineMet <- c('ALDH1A3', 'ALDH3A1', 'ALDH3B1', 'ALDH3B2', 'AOC2', 'AOC3', 'DDC', 'GOT1', 'GOT2', 'HPD', 'IL4I1', 'MAOA', 'MAOB', 'MIF', 'NAT6', 'PAH', 'PRDX6', 'TAT')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(phenylalanineMet), name = "phenylalanineMet", search=TRUE)

#33
taurineHypotaurMet <- c('ADO', 'BAAT', 'CDO1', 'CSAD', 'GAD1', 'GAD2', 'GGT1', 'GGT5', 'GGT6', 'GGT7')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(taurineHypotaurMet), name = "taurineHypotaurMet", search=TRUE)

#34
aminoacyltRNAbios <- c('AARS', 'AARS2', 'CARS', 'CARS2', 'DARS', 'DARS2', 'EARS2', 'EPRS', 'FARS2', 'FARSA', 'FARSB', 'GARS', 'HARS', 'HARS2', 'IARS', 'IARS2', 'KARS', 'LARS', 'LARS2', 'MARS', 'MARS2', 'MTFMT', 'NARS', 'NARS2', 'PARS2', 'PSTK', 'QARS', 'RARS', 'RARS2', 'SARS', 'SARS2', 'SEPSECS', 'TARS', 'TARS2', 'TARSL2', 'VARS', 'VARS2', 'WARS', 'WARS2', 'YARS', 'YARS2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(aminoacyltRNAbios), name = "aminoacyltRNAbios", search=TRUE)

#35
tryptophanMet <- c('AADAT', 'AANAT', 'ABP1', 'ACAT1', 'ACAT2', 'ACMSD', 'AFMID', 'ALDH1B1', 'ALDH2', 'ALDH3A2', 'ALDH7A1', 'ALDH9A1', 'AOX1', 'ASMT', 'CAT', 'CYP1A1', 'CYP1A2', 'CYP1B1', 'DDC', 'ECHS1', 'EHHADH', 'GCDH', 'HAAO', 'HADH', 'HADHA', 'IDO1', 'IDO2', 'IL4I1', 'INMT', 'KMO', 'KYNU', 'MAOA', 'MAOB', 'OGDH', 'OGDHL', 'TDO2', 'TPH1', 'TPH2', 'WARS', 'WARS2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(tryptophanMet), name = "tryptophanMet", search=TRUE)

#36
tyrosineMet <- c('ADH1A', 'ADH1B', 'ADH1C', 'ADH4', 'ADH5', 'ADH6', 'ADH7', 'ALDH1A3', 'ALDH3A1', 'ALDH3B1', 'ALDH3B2', 'AOC2', 'AOC3', 'AOX1', 'COMT', 'DBH', 'DCT', 'DDC', 'FAH', 'GOT1', 'GOT2', 'GSTZ1', 'HEMK1', 'HGD', 'HPD', 'IL4I1', 'LCMT1', 'LCMT2', 'MAOA', 'MAOB', 'METTL2B', 'METTL6', 'MIF', 'NAT6', 'PNMT', 'TAT', 'TH', 'TPO', 'TRMT11', 'TYR', 'TYRP1', 'WBSCR22')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(tyrosineMet), name = "tyrosineMet", search=TRUE)

#37
valineLeucineIsoleucBios <- c('BCAT1', 'BCAT2', 'IARS', 'IARS2', 'LARS', 'LARS2', 'PDHA1', 'PDHA2', 'PDHB', 'VARS', 'VARS2')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(valineLeucineIsoleucBios), name = "valineLeucineIsoleucBios", search=TRUE)

#38
Glutamine <- c('GLS2','RIMKLA','GLUD2','PYCR3','GLUL','OAT','GLUD1','GOT2','ALDH18A1','GLS','PYCR2','RIMKLB','PYCR1','KYAT1')

OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(Glutamine), name = "glutamine", search = TRUE)

#39
#Purine metabolism (denovo and salvage)
PuDeNovogenes <- c('IMPDH2','ADSS2','ATIC','PAICS','IMPDH1','GART','GMPS','LHPP','ADSS1','PFAS','PPAT','ADSL','ADSS','ADSSL1','PKM2')
OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(PuDeNovogenes), name = "PuDeNovo", search=TRUE)

#40
PuSalvagegenes <- c('AMPD1','HPRT1','APRT','GMPR2','DGUOK','DCK','PNP','GMPR','ADK','ADA','AMPD2','AMPD3','ADAL')
OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features = list(PuSalvagegenes), name = "PuSalvage", search=TRUE)

#41
#Pyrimidine metabolism (denovo and salvage)
PyDeNovogenes <- c('CAD','UMPS','DHODH')
OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features= list(PyDeNovogenes), name = "PyDeNovo", search=TRUE)

#42
PySalavagegenes <- c('UPP2','DCK','TYMP','UCK1','UCKL1','CDA','PUDP','UPP1','TK2','UCK2')
OVCA3Grpd <- AddModuleScore(OVCA3Grpd, features= list(PySalavagegenes), name= "PySalvage", search=TRUE)

saveRDS(OVCA3Grpd, paste(outdir,"OVCA3_S+P_SeuratObj-with-Metabolic-Pathway-Scores",date,".rds",sep=""))

#check for the names of the pathways
names(OVCA3Grpd[[]])

#create an object with the names of the pathways 
scores.1 <- c("OXPHOS1","glycolysis1","citricAcid1","ascorbateAldarate1","chondroitinSulfate1","fructoseMannose1","galactose1","keratanSulfate1",
              "heparanSulfate1","pentosePhosphate1","pyruvateMet1","starchSucrose1","etherLipidMet1","arachidonicAcidMet1","fattyAcidMet1","glycerolipidMet1",
              "glycerophospholipidMet1","glyoxylateDicarbMet1","linoleicAcidMet1","sphingolipidMet1","propanoateMet1","pyrimidineMet1","purineMet1",
              "alaniAspratGlutamMet1","arginineProlineMet1","butanoateMet1","glutathioneMet1","glycineSerineThreonMet1","histidineMet1","lysineDeg1",
              "cysteineMethionineMet1","phenylalanineMet1","taurineHypotaurMet1","aminoacyltRNAbios1","tryptophanMet1","tyrosineMet1",
              "valineLeucineIsoleucBios1","glutamine1","PuDeNovo1","PuSalvage1","PyDeNovo1","PySalvage1")

#export the scores to data frame
data_scores <- as.data.frame(as.matrix(OVCA3Grpd[[scores.1]]))

#save dataframe as csv
write.csv(data_scores, sep=",","~/capsule/results/metabolic_pathway_analysis/OVCA3_metabolic_pathway_scores.csv")






