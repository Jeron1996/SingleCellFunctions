## Commands to get NonNormalized, percent.mt regressed out Seurat Object for the Hansbro collaboration
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
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/master/191202-emptydrop.R")
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/master/191211-hashing_workflow.R")
source_url("https://raw.githubusercontent.com/Jeron1996/SingleCellFunctions/master/191209-QC_plots.R")
set.seed(160396)
###Start workflow by filtering out empty droplets using e.drop command
raw_data_dirs <- paste0("R://Zilog-Cancer-TumourDevelopment/190510-HANSBRO1-SingleCellData/collection", c(1:6), "/RNA/output/outs/raw_feature_bc_matrix/")
pro.name <- paste0("Collection", c(1:6))
save.Dir <- paste0("G://Hansbro_data/ReRun_normalized/Emptydrop_Out/")
save.Dir <- paste0("/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized_normalized/Emptydrop_Out/")
#Performs emptydrop for every collection
for(dir in raw_data_dirs){
  name <- strsplit(x = dir, split = "/")[[1]][5]
  raw <- read10xCounts(samples = dir, col.names = T)
  e.drop_out <- e.drop(raw_matrix_counts = raw, project.name = name, saveDir = save.Dir)
  saveRDS(object = e.drop_out, file = paste0(save.Dir, "/", name, "_e.drop_out.RDS"))
}

e.files <- list.files(path = "/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized_normalized/Emptydrop_Out/" )[grepl(x = list.files(path = "/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized/Emptydrop_Out/" ), pattern = "_e.drop_out.RDS", ignore.case = T)]

#Translation for hashtags, need for hashing_workflow
hash <- c("A0303", "A0304", "A0305", "A0306", "A0307", "A0308")
hash_trans <- as.character(c("Sal/SPG", "Ova 1", "Ova/DEX", 'Sal/Cmu', "Ova/Cmu", "Ova/Cmu/DEX"))
hashing_translation <- data.frame(hash, hash_trans, stringsAsFactors = F)
output_dir <- "/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized/HashingOut/"

#Loops over all files of the emptydrop output and performs the hashing workflow to translate hashtags into samples and filter out droplets with double hashtags.
for(file in e.files){
  e.out <- readRDS(file = paste0(save.Dir, file))
  projectName <- strsplit(x = file, split = "_")[[1]][1]
  umi_directory <- paste0("/share/ScratchGeneral/jerven/Hansbro_data/Data/", projectName, "_output/umi_count/")
  hto_directory <- paste0("/share/ScratchGeneral/jerven/Hansbro_data/Data/", projectName, "_output/read_count/")
  hashtag_output <- cell_hashing_workflow(seurat.obj.emptydrop = e.out, plotting = TRUE, saveDir = output_dir, umi_dir = umi_directory, hto_dir = hto_directory, pro.name = projectName, hash_translation = hashing_translation)
  hashtag_output@project.name <- projectName
  saveRDS(object = hashtag_output, file = paste0(output_dir, "/", projectName, "_hashingOutWithSeed.RDS"))
}

seurat_list <- list()
hashing_files <- list.files("/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized/HashingOut/")[grepl(x = list.files("/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized/HashingOut/"), pattern = "WithSeed.RDS")]

#create list with all hashing output files, to be used for merging
for(hashing_file in hashing_files){
  cell_prefix <- strsplit(x = hashing_file, split = "_")[[1]][1]
  cell_prefix <- strsplit(x = cell_prefix, split = "collection")[[1]][2]
  seurat.obj.emptydrop <- readRDS(file = paste0(output_dir, "/", hashing_file))
  seurat.obj.emptydrop <- RenameCells(object = seurat.obj.emptydrop, add.cell.id = paste0("c", cell_prefix))
  seurat_list <- append(value=seurat.obj.emptydrop, x=seurat_list)
}

#Merge all hashtag outputs
seurat_merged <- merge(x = seurat_list[[1]], y = c(seurat_list[[2]], seurat_list[[3]], seurat_list[[4]], seurat_list[[5]], seurat_list[[6]]))
seurat_merged <- seurat_merged[, seurat_merged$percent.mt < 10]
seurat_dir <- "/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized/SeuratObjects/"
plot_dir <- "/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized/Plots/"

#Perform simple Seurat workflow, make plots that can be used for reference.
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged)
seurat_merged <- ScaleData(object = seurat_merged, features = rownames(seurat_merged))
seurat_merged <- RunPCA(seurat_merged, features = VariableFeatures(object = seurat_merged))
#Use first 20 dimensions for UMAP
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:15)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:15)
#Add Cell cycle information to the Seurat Object
load("/share/ScratchGeneral/jerven/Hansbro_data/CellCylceGenes/cell.cyclegenes.rdata")
seurat_merged <- CellCycleScoring(object = seurat_merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

#Save Seurat Object
saveRDS(object = seurat_merged, file = paste0(seurat_dir, "/200227_ALL_merged_normalized.RDS"))

#Make and save several Plots
plots_cluster(seurat.object = seurat_merged, save.name = "200227_ALL_merged_normalized", dir = plot_dir)

#Perform Seurat analysis again, but this time regress out percent.mt

seurat_regressed <- readRDS(file = paste0(seurat_dir, "/200227_ALL_merged_normalized.RDS"))
resolution <- c(0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.2, 1.4)
seurat_regressed <- NormalizeData(seurat_regressed)
seurat_regressed <- ScaleData(seurat_regressed, vars.to.regress = "percent.mt", features = rownames(seurat_regressed))
#Add Cell cycle information to the Seurat Object
load("/share/ScratchGeneral/jerven/Hansbro_data/CellCylceGenes/cell.cyclegenes.rdata")
seurat_regressed <- CellCycleScoring(object = seurat_regressed, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
seurat_regressed <- FindVariableFeatures(seurat_regressed)
VariableFeatPlot <- VariableFeaturePlot(seurat_regressed)
ggsave(plot=VariableFeatPlot, filename="/share/ScratchGeneral/jerven/Hansbro_data/ReRun_normalized/plots/200227_ALL_merged_MT_regressed_VariableFeaturePlot.pdf", width=10, height=10)
seurat_regressed <- RunPCA(seurat_regressed, features = VariableFeatures(object = seurat_regressed))
seurat_regressed <- FindNeighbors(seurat_regressed, dims = 1:15)
for(res in resolution){
  seurat_regressed <- FindClusters(seurat_regressed, resolution = res)
}
seurat_regressed <- RunUMAP(seurat_regressed, dims = 1:15)

#Save Regressed seurat output
saveRDS(object = seurat_regressed, file = paste0(seurat_dir, "/200227_ALL_merged_MT_regressed_normalized.RDS"))

#Make plots
plots_cluster(seurat.object = seurat_regressed, save.name = "200227_ALL_merged_MT_regressed_normalized", dir = plot_dir)
