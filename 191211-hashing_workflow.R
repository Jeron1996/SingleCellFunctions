### Simple cell hashing workflow to identify the number of hash-tags and the corresponding number of cells

library(Seurat)
library(ggplot2)
library(dplyr)
library(Matrix)
library(stringr)

cell_hashing_workflow <- function(plotting = TRUE, pos.quart = 0.99, seurat.obj.emptydrop, saveDir, umi_dir, pro.name, hash_translation){
  set.seed(160396)
  ##Creates UMI SeuratObject from UMI count matrix and formats it
  UMI <- Read10X(data.dir = umi_dir, gene.column = 1)
  rownames(UMI) <- substr(x = rownames(UMI), start = 1, stop = 5)
  UMI <- CreateSeuratObject(counts = UMI, project = pro.name, assay = "HTO")
  UMI <- UMI[rownames(UMI)[rownames(UMI) != "unmap"],]

  ##filter empty droplets based on e.drops output
  ok.barcodes <- colnames(seurat.obj.emptydrop)
  ok.barcodes.clean <- str_remove(string = ok.barcodes, pattern = "-1")
  UMI_empty <- UMI[, ok.barcodes.clean]
  n.ok.cells.edrop <- length(ok.barcodes.clean)
  n.ok.cells.hash <- length(colnames(UMI_empty))
  print(paste0("Number of ok cells: ", n.ok.cells.edrop, "; Number of cells in UMI: ", n.ok.cells.hash))

  ##Seurat hashing workflow
  hashtag <- UMI_empty

  ##Normalize using log-normalization
  hashtag <- NormalizeData(object = hashtag, assay = "HTO", normalization.method = "CLR", verbose = FALSE)

  #optional using SCTransform normalization instead of log-normalization
  #hashtag <- SCTransform(hashtag)

  ##Identify predominat hash.tag for every cell
  ###play with init, pos.quantile and k.func parameters in HTODemux to tweak selection cutoffs for label calling
  ###Doesnt matter if we use SCTransform normalization or NormalizeData, did not affect the hash.tag counts in test dataset
  hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = pos.quart)

  ##Print number of cells annotated as singlet (one hash.tag), doublets (contain multiple hash.tags) or negative (no hash.tag)
  t <- table(hashtag$HTO_classification.global)
  print(t)

  ##Preform hashing without emptydrop filtering
  hashtag_all <- UMI
  hashtag_all <- NormalizeData(object = hashtag_all, assay = "HTO", normalization.method = "CLR", verbose = FALSE)
  hashtag_all <- HTODemux(hashtag_all, assay = "HTO", positive.quantile = pos.quart)
  saveRDS(object=hashtag_all, file = paste0(hashing.dir, "/", projectName, "_hashingNoEmpty.RDS"))
  
  ##create QC plots, if desired
  if(plotting){
    saveDir <- saveDir
    QC_plots_hashing(hashtag = hashtag, pos.quart = pos.quart, saveDir = saveDir, name = pro.name)
    }

  ##Add hastag annotation to emptydrop output
  colna_hash <- paste0(colnames(hashtag), "-1")
  colna_hash <- sort(colna_hash)
  colna_empty <- colnames(seurat.obj.emptydrop)
  colna_empty <- colna_empty[is.element(colna_empty, colna_hash)]
  colna_empty <- sort(colna_empty)
  ord.log <- colna_hash == colna_empty

  if(sum(ord.log) == length(colna_empty)){
    seurat.obj.emptydrop <- seurat.obj.emptydrop[, colna_hash]
    } else{
      stop("Barcodes not identical or in the same order")
    }

  ##Add hashtag data to emptydrop output
  hashtag <- RenameCells(object = hashtag, new.names = colna_hash)
  seurat.obj.emptydrop <- AddMetaData(object = seurat.obj.emptydrop, metadata = hashtag@meta.data, col.name = c(colnames(hashtag@meta.data)))

  ##Return only singlets & translate hashtags
  seurat.obj.emptydrop <- subset(x = seurat.obj.emptydrop, subset = HTO_classification.global %in% "Singlet")
  seurat.obj.emptydrop$SampleID <- NA
  for(i in hash_translation$hash){
    seurat.obj.emptydrop$SampleID[seurat.obj.emptydrop$HTO_classification %in% i] <- hash_translation[hash_translation$hash %in% i, "hash_trans"]
  }

  seurat.obj.emptydrop
}
