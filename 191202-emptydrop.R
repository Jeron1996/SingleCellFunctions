#######EmptyDrop workflow

library(DropletUtils)
library(Seurat)
library(ggplot2)


#optimze so that all p-values at/below limit are uniformly distributed


e.drop <- function(raw_matrix_counts, project.name, saveDir){
  set.seed(160396)
  data <- raw_matrix_counts
  gene.symbols <- rowData(data)[2]
  gene.symbols <- gene.symbols$Symbol
  rownames(data) <- gene.symbols

  #Convert raw data to counts
  mycounts <- counts(data)

  ##filter out genes which are not found in any barcode
  detected <- rowSums(x = mycounts) == 0

  if(sum(is.na(detected)) != 0){
    stop("NAs detected when filtering genes")
  }

  mycounts.detect <- mycounts[!detected, , drop = F]

  ##Make barcodePlot
  #barcodeplot(counts = mycounts.detect)

  ## Run emptyDrops()
  print("This is going to take a while")
  e.out <- emptyDrops(m = mycounts.detect, lower = 100, ignore = 1)
  print("This took a while")
  print(e.out)

  ##Create dgcMatrix with cell and empty droplets

  is.cell <- e.out$FDR <= 0.05
  print("if sig == FALSE and limited == TRUE is > 0, increase niter in emptydrops")
  print(table(Limited=e.out$Limited, Significant=is.cell))

  ok.cells <- sum(is.cell, na.rm=TRUE)                     #Shows the number of cells
  print(paste0("Ok cells in data: ", ok.cells))

  is.no.cell <- !is.cell
  empty.cells <- sum(is.no.cell, na.rm=TRUE)
  print(paste0("Empty cells in data: ", empty.cells))

  qc_value <- sum(is.cell == is.no.cell, na.rm = T)  ##should be 0, indicates that no values for is.cell and is.no.cell are the same
  print(paste0("Difference between is.cell and is.no.cell: ", qc_value))

  is.cell[is.na(is.cell)] <- FALSE
  is.no.cell[is.na(is.no.cell)] <- FALSE

  data.empty <- data[, is.no.cell]
  data.ok <- data[, is.cell]

  counts.ok <- counts(data.ok)
  counts.empty <- counts(data.empty)

  ##Create Seurat objects from new data matrix files
  seurat.ok <- CreateSeuratObject(counts = counts.ok, project = as.character(project.name), min.cells = 3)
  seurat.empty <- CreateSeuratObject(counts = counts.empty, project = as.character(project.name), min.cells = 3)

  #Calculate percentage mitochondrial genes
  seurat.ok[["percent.mt"]] <- PercentageFeatureSet(seurat.ok, pattern = "^mt-")
  seurat.ok <- subset(seurat.ok, subset = nFeature_RNA <= 4000)   #Cut off all cells with more than 4000 features
  seurat.empty[["percent.mt"]] <- PercentageFeatureSet(seurat.empty, pattern = "^mt-")

  saveDir <- saveDir

  QC_plots_emptydrop(seuratObj.ok = seurat.ok, seuratObj.empty = seurat.empty, pro.name = project.name, low = 100, FDR = 0.05, saveDir = saveDir, eOut = e.out)

  seurat.ok
}
