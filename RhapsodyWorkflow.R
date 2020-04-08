## Commands to get SCT normalized, percent.mt regressed out Seurat Object for the Hansbro collaboration
library(ggforce)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(stringr)
library(DropletUtils)
library(dplyr)
library(Matrix)
library(devtools)
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/Rhapsody/QC_plots_Rhapsody.R")
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/Rhapsody/hashing_rhapsody.R")
set.seed(160396)

##Choose normalization method
norm.methods <- c("SCTransform", "log-transform", "none")
norm.method <- norm.methods[1]
print(paste0("Normalization method: ", norm.method))

##Create a new output directory for this analysis run
date <- Sys.Date()
date <- format(date, "%y%m%d")
ProjectTitle <- "Nova_bothLanes"
ProjectName <- paste0(date, "_", ProjectTitle)
if(dir.exists(paste0("/share/ScratchGeneral/jerven/Rhapsody/Analysis/", ProjectName))){
  stop("Directory belonging to this ProjectName already exists. \n  Please choose another ProjectName")
}
dir.path <- paste0("/share/ScratchGeneral/jerven/Rhapsody/Analysis", ProjectName, "/", c("Plots", "SeuratObjects", "DEgenes"))
sapply(dir.path, function(X) dir.create(path = X, recursive = T))

##Directories and variables needed
raw_data_dirs <- "/share/ScratchGeneral/jerven/Rhapsody/Nova/BothLanes/_1_Combined_AAGAGGCA_DBEC_MolsPerCell.csv"
pro.name <- "TestRun"
plot_dir <- dir.path[1]
seurat_dir <- dir.path[2]
de_dir <- dir.path[3]

##Load UMIs per cell file & format so that it can be used as Seurat input
umi_reads <- read.csv(file = raw_data_dirs, comment.char="#", stringsAsFactors=FALSE)
rownames(umi_reads) <- umi_reads[,1]
umi_reads <- umi_reads[, -1]          #Set cell identeties as rownames and delete that column from dataset
gene_merged <- strsplit(colnames(umi_reads), split="_") #Gene IDs are merged together into one string by SB Rhaspsody. This splits them
gene_merged[sapply(gene_merged, FUN = function(x) x[2] == "BD")][[1]][1] <- paste0(gene_merged[sapply(gene_merged, FUN = function(x) x[2] == "BD")][[1]][1], "-BD") #Some genes are BD specific, change their name
GeneSymbol <- sapply(X=gene_merged, FUN = function(X) X[1]) #Select GeneSymbols
ENSMBL_ID <- sapply(X=gene_merged, FUN = function(X) X[2]) #Select Ensmlb IDs
colnames(umi_reads) <- GeneSymbol #Set Gene symbols as gene names
umi_reads <- t(umi_reads) #transpose counts, so that it can be read by Seurat

umi_Reads <- CreateSeuratObject(counts = umi_reads, project = ProjectTitle, min.cells = 1, min.features = 1)

#Add sample demultiplexing information
umi_Reads <- hashing_rhapsody(rhap_obj = umi_Reads)

##Choose which normalization method to perform on the merged data.
if(norm.method == "SCTransform"){ #Perform SCTransform normalization
  umi_Reads <- SCTransform(object = umi_Reads, verbose = TRUE)
}else if(norm.method == "log-transform"){ #Perform log normalization
  umi_Reads <- NormalizeData(umi_Reads)
  umi_Reads <- FindVariableFeatures(umi_Reads)
  umi_Reads <- ScaleData(object = umi_Reads, features = rownames(umi_Reads))
}else if(norm.method == "none"){ #Perform NO normalization, not advised
  umi_Reads <- FindVariableFeatures(umi_Reads)
  umi_Reads <- ScaleData(object = umi_Reads, features = rownames(umi_Reads))
}else{
  print("No normalization method chosen")
}

#RunPCA using standard values and do UMAP
umi_Reads <- RunPCA(umi_Reads)
umi_Reads <- RunUMAP(umi_Reads, dims = 1:25)

#Use first 25 dimensions for UMAP
umi_Reads <- FindNeighbors(umi_Reads, dims = 1:25)

resolution <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.2, 1.4)
for(res in resolution){
  umi_Reads <- FindClusters(umi_Reads, resolution = res)
}

#Save Seurat Object
saveRDS(object = umi_Reads, file = paste0(seurat_dir, "/", ProjectName, ".RDS"))

#save objects for SingleR workflow on local machine
expression <- GetAssayData(object = umi_Reads, assay = "SCT", slot = "data")
saveRDS(object = expression, file = paste0(seurat_dir, "/SCT_data_expression.RDS"))

#Make and save several Plots
plots_cluster(seurat.object = umi_Reads, save.name = ProjectName, dir = plot_dir)

##Find Marker genes associated to all resolutions
resol <- umi_Reads@meta.data
resol <- resol[grepl(pattern = "res.", x=resol)]
for(re in resol){
  Idents(umi_Reads) <- re
  markers <- FindAllMarkers(object = umi_Reads, assay = "SCT", slot = "data", logfc.threshold = 0.25)
  saveRDS(object = markers, file=paste0(de_dir, "/", re, "_markers.RDS"))
}
