###Quality Control Summary
###Creates multiple plots to show the quality of the emptydrop result

library(ggplot2)
library(Seurat)
library(clustree)

###################################################################################################################################
####################################################QC plots emptydrop workflow####################################################
###################################################################################################################################

QC_plots_emptydrop <- function(seuratObj.ok, seuratObj.empty, low = low, pro.name, FDR, saveDir){

  ##Format Seurat Objects & Calculate percent mt
  seuratObj.ok$orig.ident <- "Ok"
  seuratObj.empty$orig.ident <- "Empty"

  seurat.combine <- merge(x = seuratObj.ok, y = seuratObj.empty, add.cell.ids = c("Ok", "Empty"), project = "QC_plots")

  ##Create boxplot for nCount_RNA
  df.count.ok <- data.frame(seuratObj.ok$nCount_RNA, "Ok")
  df.count.empty <- data.frame(seuratObj.empty$nCount_RNA, "Empty")

  colnames(df.count.ok) <- c("Counts", "Droplet_state")
  colnames(df.count.empty) <- c("Counts", "Droplet_state")

  df.count.comb <- rbind(df.count.ok, df.count.empty)
  df.count.comb$Droplet_state <- factor(df.count.comb$Droplet_state, levels = c("Empty", "Ok"))
  xlabel.count <- paste0("Number of ok cells: ", nrow(df.count.comb[df.count.comb$Droplet_state == "Ok", ]),
                         "; Number of empty droplets: ", nrow(df.count.comb[df.count.comb$Droplet_state == "Empty",]))
  box_count <- ggplot(data = df.count.comb, aes(x = Droplet_state, y = Counts)) + geom_boxplot() + facet_zoom(ylim = c(0, 500)) +
                          theme_bw() + xlab(label = xlabel.count)

  ##Create boxplot for nFeature_RNA
  df.feat.ok <- data.frame(seuratObj.ok$nFeature_RNA, "Ok")
  df.feat.empty <- data.frame(seuratObj.empty$nFeature_RNA, "Empty")

  colnames(df.feat.ok) <- c("Features", "Droplet_state")
  colnames(df.feat.empty) <- c("Features", "Droplet_state")

  df.feat.comb <- rbind(df.feat.ok, df.feat.empty)
  df.feat.comb$Droplet_state <- factor(df.feat.comb$Droplet_state, levels = c("Empty", "Ok"))

  xlabel.feat <- paste0("Number of ok cells: ", nrow(df.feat.comb[df.feat.comb$Droplet_state == "Ok", ]),
                        "; Number of empty droplets: ", nrow(df.feat.comb[df.feat.comb$Droplet_state == "Empty",]))
  box_feat <- ggplot(data = df.feat.comb, aes(x = Droplet_state, y = Features)) + geom_boxplot() +
                        facet_zoom(ylim = c(0, 500)) + theme_bw() +   xlab(label = xlabel.feat)


  ### Make ViolinPlot comparing mitochondrial gene percentage
  Vln_mt <- VlnPlot(object = seurat.combine, features = "percent.mt", split.by = "orig.ident", pt.size = 0)


  ### Save plots in single PDF file
  pdf.name <- paste0(saveDir, "/", pro.name,"_lower-", low, "_FDR-", FDR, "_emptydrop_QC.pdf")
  pdf(file = pdf.name)
  plot(box_feat)
  plot(box_count)
  plot(Vln_mt)
  dev.off()

  boxplot.list <- list(df.feat.comb, df.feat.comb)
  names(boxplot.list) <- c("boxplot_data_features", "boxplot_data_counts")
  list.name <- paste0(saveDir, pro.name,"_lower-", low, "_FDR-", FDR, "_boxplotData.rds")
  saveRDS(object = boxplot.list, file = list.name)
  print("Plots saved")
}


###################################################################################################################################
#####################################################QC plots hashing workflow#####################################################
###################################################################################################################################

library(Seurat)
library(ggplot2)
library(gridExtra)
library(stringr)

##Make contigency table, ridge plot and tSNE plot to show counts per hash tag
QC_plots_hashing <- function(hashtag = hashtag, pos.quart = pos.quart, saveDir, name){

  #Idents(hashtag) <- "HTO_maxID"

  ##make ridge plot for to show expression of each hashtag
  ridge.plot <- RidgePlot(hashtag, assay = "RNA", features = rownames(hashtag[["RNA"]]), ncol = 2)

  ##make heatmap of expression values
  heatmap.hashtag <- HTOHeatmap(hashtag, assay = "HTO", ncells = 5000)

  ## make tSNE plot for all hashtags except for negative cells
  hashtag.subset <- subset(hashtag, idents = "Negative", invert = TRUE)
  #rownames(hashtag.subset) <- substr(x = rownames(hashtag.subset), start = 1, stop = 5)
  #hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = hashtag.subset, assay = "HTO"))))
  #hashtag.subset <- RunTSNE(hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
  #dim.plot <- DimPlot(hashtag.subset, label = T, point.size = 1)
  #dim.plot.global <- DimPlot(hashtag.subset, label = T, group.by = "HTO_classification.global", point.size = 1)

  ##format contingency table
  count_table <- table(hashtag$HTO_classification)
  count_df <- as.data.frame(count_table)
  count_df <- count_df[base::order(x = count_df[, 2], decreasing = TRUE), ]
  count_df <- data.frame(count_df[1:(nrow(count_df)/2), ], count_df[((nrow(count_df)/2)+1):nrow(count_df), ])
  colnames(count_df) <-c("Hashtag" ,"Frequency", "Hashtag" ,"Frequency")

  ##save as PDF
  project.name <- hashtag@project.name
  project.name <- str_replace(project.name, "[[:punct:]]", ",")
  pos.quart <- str_replace(pos.quart, "[[:punct:]]", ",")

  ##save as pdf
  pdf(file = paste0(saveDir, "/", name, "_QuartileCutOff-", pos.quart, "_Plots.pdf"), paper="a4")
  grid.table(count_df)
  plot(heatmap.hashtag)
  plot(ridge.plot)
  dev.off()
}

###################################################################################################################################
########################################################QC plots clustering########################################################
###################################################################################################################################

plots_cluster <- function(seurat.object, save.name, dir){
  obj <- seurat.object
  plot_dir <- dir
  sv.name <- save.name

  ## Generate and save plots
  Nelbow <- ElbowPlot(object = obj, ndims = 30)
  ggsave(filename = paste0(plot_dir, "/", sv.name, "_withSeed_ElbowPlot.pdf"), plot = Nelbow, device = "pdf", width = 10, height = 10)

  NsampleID <- DimPlot(object = obj, group.by = "SampleID", label = T)
  ggsave(filename = paste0(plot_dir, "/", sv.name, "_SampleID.pdf"), plot = NsampleID, device = "pdf", width = 10, height = 10)

  Norig.ident <- DimPlot(object = obj, group.by = "orig.ident", label = T)
  ggsave(filename = paste0(plot_dir, "/", sv.name, "_OrigIdent.pdf"), plot = Norig.ident, device = "pdf", width = 10, height = 10)

  NCellCycle <- DimPlot(object = obj, group.by = "Phase", label = T)
  ggsave(filename = paste0(plot_dir, "/", sv.name, "_CellCyclePhase.pdf"), plot = NCellCycle, device = "pdf", width = 10, height = 10)

  NpercentMT <- FeaturePlot(object = obj, features = "percent.mt")
  ggsave(filename = paste0(plot_dir, "/", sv.name, "_PercentMT.pdf"), plot = NpercentMT, device = "pdf", width = 10, height = 10)

  ## Plot DimPlots for all resolutions
  resolution <- names(obj@meta.data)
  resolution <- resolution[grepl(pattern = "res.", x = resolution)]

  if(length(resolution) != 0){
    for(res in resolution){
      Idents(obj) <- res

      temp_plot <- DimPlot(obj, reduction="umap", label=T)
      ggsave(filename = paste0(plot_dir, "/", sv.name, "_", res, ".pdf"), plot = temp_plot, device = "pdf", width = 10, height = 10)

      nCount.vln <- VlnPlot(object=obj, features="nCount_RNA", pt.size=0.001) + ggtitle("nCount_RNA distribution per cluster")
      ggsave(filename = paste0(plot_dir, "/", sv.name, "_", res, "_nCount.pdf"), plot = nCount.vln, device = "pdf", width = 10, height = 10)

      nFeature.vln <- VlnPlot(object=obj, features="nFeature_RNA", pt.size=0.001) + ggtitle("nFeature_RNA distribution per cluster")
      ggsave(filename = paste0(plot_dir, "/", sv.name, "_", res, "_nFeature.pdf"), plot = nFeature.vln, device = "pdf", width = 10, height = 10)

      percent.mt.vln <- VlnPlot(object=obj, features="percent.mt", pt.size=0.001) + ggtitle("Percentage mitochondrial genes per cluster")
      ggsave(filename = paste0(plot_dir, "/", sv.name, "_", res, "_PercentMT_VLN.pdf"), plot =   percent.mt.vln, device = "pdf", width = 10, height = 10)
    }

    clustre <- clustree(x=seurat, prefix="RNA_snn_res.")
    ggsave(filename = paste0(plot_dir, "/", sv.name, "_Clustree.pdf"), plot = clustre, device = "pdf", width = 10, height = 10)
  }

}
