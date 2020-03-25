## Commands to get SCT normalized, percent.mt regressed out Seurat Object for the Hansbro collaboration
.libPaths("/share/ClusterShare/software/contrib/jerven/R/3.6/")
library(ggforce)
library(Seurat)
library(ggplot2)
library(gridExtra)
library(stringr)
library(DropletUtils)
library(dplyr)
library(Matrix)
library(devtools)
<<<<<<< Updated upstream
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/master/191202-emptydrop.R")
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/master/191211-hashing_workflow.R")
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/master/191209-QC_plots.R")
=======
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/swithMergeHash/191202-emptydrop.R")
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/swithMergeHash/191211-hashing_workflow.R")
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/swithMergeHash/191209-QC_plots.R")
>>>>>>> Stashed changes
set.seed(160396)

##Choose normalization method
norm.methods <- c("SCTransform", "log-transform", "none")
norm.method <- norm.methods[1]
print(paste0("Normalization method: ", norm.method))

##Create a new output directory for this analysis run
date <- Sys.Date()
date <- format(date, "%y%m%d")
ProjectTitle <- "SCT_retain_ALL"
ProjectName <- paste0(date, "_", ProjectTitle)
if(dir.exists(paste0("/share/ScratchGeneral/jerven/Hansbro_data/Analysis/", ProjectName))){
  stop("Directory belonging to this ProjectName already exists. \n  Please choose another ProjectName")
}
dir.path <- paste0("/share/ScratchGeneral/jerven/Hansbro_data/Analysis/", ProjectName, "/", c("Emptydrop_Out", "HashingOut", "Plots", "SeuratObjects", "DEgenes"))
sapply(dir.path, function(X) dir.create(path = X, recursive = T))

##Directories and variables needed
raw_data_dirs <- paste0("/share/ScratchGeneral/jerven/Hansbro_data/Data/RNA/collection", c(1:6), "/raw_feature_bc_matrix")
pro.name <- paste0("Collection", c(1:6))
emptydrop.dir <- dir.path[1]
hashing.dir <- dir.path[2]
plot_dir <- dir.path[3]
seurat_dir <- dir.path[4]
de_dir <- dir.path[5]

###Start workflow by filtering out empty droplets using e.drop command
#Performs emptydrop for every collection
for(dir in raw_data_dirs){
  name <- strsplit(x = dir, split = "/")[[1]][8]
  raw <- read10xCounts(samples = dir, col.names = T)
  e.drop_out <- e.drop(raw_matrix_counts = raw, project.name = name, saveDir = emptydrop.dir)
  saveRDS(object = e.drop_out, file = paste0(emptydrop.dir, "/", name, "_e.drop_out.RDS"))
}


e.files <- list.files(path = emptydrop.dir)[grepl(x = list.files(path = emptydrop.dir), pattern = "_e.drop_out.RDS", ignore.case = T)]

#Translation for hashtags, need for hashing_workflow
hash <- c("A0303", "A0304", "A0305", "A0306", "A0307", "A0308")
hash_trans <- as.character(c("Sal/SPG", "Ova 1", "Ova/DEX", 'Sal/Cmu', "Ova/Cmu", "Ova/Cmu/DEX"))
hashing_translation <- data.frame(hash, hash_trans, stringsAsFactors = F)

#Loops over all files of the emptydrop output and performs the hashing workflow to translate hashtags into samples and filter out droplets with double hashtags.
for(file in e.files){
  e.out <- readRDS(file = paste0(emptydrop.dir, "/", file))
  projectName <- strsplit(x = file, split = "_")[[1]][1]
  umi_directory <- paste0("/share/ScratchGeneral/jerven/Hansbro_data/Data/", projectName, "_output/umi_count/")
  hashtag_output <- cell_hashing_workflow(seurat.obj.emptydrop = e.out, plotting = TRUE, saveDir = hashing.dir, umi_dir = umi_directory, pro.name = projectName, hash_translation = hashing_translation)
  hashtag_output@project.name <- projectName
  saveRDS(object = hashtag_output, file = paste0(hashing.dir, "/", projectName, "_hashingOut.RDS"))
}

seurat_list <- list()
hashing_files <- list.files(hashing.dir)[grepl(x = list.files(hashing.dir), pattern = "_hashingOut.RDS")]

#create list with all hashing output files, to be used for merging
for(hashing_file in hashing_files){
  cell_prefix <- strsplit(x = hashing_file, split = "_")[[1]][1]
  cell_prefix <- strsplit(x = cell_prefix, split = "collection")[[1]][2]
  seurat.obj.emptydrop <- readRDS(file = paste0(hashing.dir, "/", hashing_file))
  seurat.obj.emptydrop <- RenameCells(object = seurat.obj.emptydrop, add.cell.id = paste0("c", cell_prefix))
  seurat_list <- append(value=seurat.obj.emptydrop, x=seurat_list)
}





e.files <- list.files(path = emptydrop.dir)[grepl(x = list.files(path = emptydrop.dir), pattern = "_e.drop_out.RDS", ignore.case = T)]

#Translation for hashtags, need for hashing_workflow
hash <- c("A0303", "A0304", "A0305", "A0306", "A0307", "A0308")
hash_trans <- as.character(c("Sal/SPG", "Ova 1", "Ova/DEX", 'Sal/Cmu', "Ova/Cmu", "Ova/Cmu/DEX"))
hashing_translation <- data.frame(hash, hash_trans, stringsAsFactors = F)

#Loops over all files of the emptydrop output and performs the hashing workflow to translate hashtags into samples and filter out droplets with double hashtags.
for(file in e.files){
  e.out <- readRDS(file = paste0(emptydrop.dir, "/", file))
  projectName <- strsplit(x = file, split = "_")[[1]][1]
  umi_directory <- paste0("/share/ScratchGeneral/jerven/Hansbro_data/Data/", projectName, "_output/umi_count/")
  hashtag_output <- cell_hashing_workflow(seurat.obj.emptydrop = e.out, plotting = TRUE, saveDir = hashing.dir, umi_dir = umi_directory, pro.name = projectName, hash_translation = hashing_translation)
  hashtag_output@project.name <- projectName
  saveRDS(object = hashtag_output, file = paste0(hashing.dir, "/", projectName, "_hashingOut.RDS"))
}

seurat_list <- list()
empty_files <- list.files(emptydrop.dir)[grepl(x = list.files(emptydrop.dir), pattern = "_e.drop_out.RDS")]

#create list with all hashing output files, to be used for merging
for(empty_file in empty_files){
  cell_prefix <- strsplit(x = empty_file, split = "_")[[1]][1]
  cell_prefix <- strsplit(x = cell_prefix, split = "collection")[[1]][2]
  seurat.obj.emptydrop <- readRDS(file = paste0(emptydrop.dir, "/", empty_file))
  seurat.obj.emptydrop <- RenameCells(object = seurat.obj.emptydrop, add.cell.id = paste0("c", cell_prefix))
  seurat_list <- append(value=seurat.obj.emptydrop, x=seurat_list)
}






#Merge all hashtag outputs
seurat_merged <- merge(x = seurat_list[[1]], y = c(seurat_list[[2]], seurat_list[[3]], seurat_list[[4]], seurat_list[[5]], seurat_list[[6]]))
seurat_merged <- seurat_merged[, seurat_merged$percent.mt < 10]

##Choose which normalization method to perform on the merged data.
if(norm.method == "SCTransform"){ #Perform SCTransform normalization
  seurat_merged <- SCTransform(object = seurat_merged, verbose = TRUE)
}else if(norm.method == "log-transform"){ #Perform log normalization
  seurat_merged <- NormalizeData(seurat_merged)
  seurat_merged <- FindVariableFeatures(seurat_merged)
  seurat_merged <- ScaleData(object = seurat_merged, features = rownames(seurat_merged))
}else if(norm.method == "none"){ #Perform NO normalization, not advised
  seurat_merged <- FindVariableFeatures(seurat_merged)
  seurat_merged <- ScaleData(object = seurat_merged, features = rownames(seurat_merged))
}else{
  print("No normalization method chosen")
}

#RunPCA using standard values
seurat_merged <- RunPCA(seurat_merged)
#Use first 25 dimensions for UMAP
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:25)

resolution <- c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.2, 1.4)
for(res in resolution){
  seurat_merged <- FindClusters(seurat_merged, resolution = res)
}

seurat_merged <- RunUMAP(seurat_merged, dims = 1:25)
#Add Cell cycle information to the Seurat Object
load("/share/ScratchGeneral/jerven/Hansbro_data/CellCylceGenes/cell.cyclegenes.rdata")
seurat_merged <- CellCycleScoring(object = seurat_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

#Save Seurat Object
saveRDS(object = seurat_merged, file = paste0(seurat_dir, "/", ProjectName, ".RDS"))

#Make and save several Plots
plots_cluster(seurat.object = seurat_merged, save.name = ProjectName, dir = plot_dir)

##Find Marker genes associated to all resolutions
resol <- seurat_merged@meta.data
resol <- resol[grepl(pattern = "res.", x=resol)]
for(re in resol){
  Idents(seurat_merged) <- re
  markers <- FindAllMarkers(object = seurat_merged, assay = "RNA", logfc.threshold = 0.25)
  saveRDS(object = markers, file=paste0(de_dir, "/", re, "_markers.RDS"))
}

#Perform Seurat analysis again, but this time regress out percent.mt

seurat_regressed <- readRDS(file = paste0(seurat_dir, "/", ProjectName, ".RDS"))

##Choose which normalization method to perform on the merged data.
if(norm.method == "SCTransform"){ #Perform SCTransform normalization
  seurat_regressed <- SCTransform(object = seurat_regressed, vars.to.regress = "percent.mt", verbose = TRUE)
}else if(norm.method == "log-transform"){ #Perform log normalization
  seurat_regressed <- NormalizeData(seurat_regressed)
  seurat_regressed <- FindVariableFeatures(seurat_regressed)
  seurat_regressed <- ScaleData(object = seurat_regressed, features = rownames(seurat_regressed))
}else if(norm.method == "none"){ #Perform NO normalization, not advised
  seurat_regressed <- FindVariableFeatures(seurat_regressed)
  seurat_regressed <- ScaleData(object = seurat_regressed, features = rownames(seurat_regressed))
}else{
  print("No normalization method chosen")
}

#Add Cell cycle information to the Seurat Object
load("/share/ScratchGeneral/jerven/Hansbro_data/CellCylceGenes/cell.cyclegenes.rdata")
seurat_regressed <- CellCycleScoring(object = seurat_regressed, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
seurat_regressed <- RunPCA(seurat_regressed)
seurat_regressed <- FindNeighbors(seurat_regressed, dims = 1:25)
for(res in resolution){
  seurat_regressed <- FindClusters(seurat_regressed, resolution = res)
}
seurat_regressed <- RunUMAP(seurat_regressed, dims = 1:25)

#Save Regressed seurat output
saveRDS(object = seurat_regressed, file = paste0(seurat_dir, "/", ProjectName, "_MT_regressed.RDS"))

#Make plots
plots_cluster(seurat.object = seurat_regressed, save.name = paste0(ProjectName, "_MT_regressed"), dir = plot_dir)

##Find Marker genes associated to all resolutions after regression
resol_regressed <- seurat_regressed@meta.data
resol_regressed <- resol_regressed[grepl(pattern = "res.", x=resol_regressed)]
for(re_regressed in resol_regressed){
  Idents(seurat_regressed) <- re_regressed
  markers <- FindAllMarkers(object = seurat_regressed, assay = "RNA", logfc.threshold = 0.25)
  saveRDS(object = markers, file=paste0(de_dir, "/", re_regressed, "MT_regressed_markers.RDS"))
}
