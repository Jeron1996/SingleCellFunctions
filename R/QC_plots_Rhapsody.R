###Quality Control Summary
###Creates multiple plots to show the quality of the emptydrop result

library(ggplot2)
library(Seurat)
library(clustree)

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

  Norig.ident <- DimPlot(object = obj, group.by = "orig.ident", label = T)
  ggsave(filename = paste0(plot_dir, "/", sv.name, "_OrigIdent.pdf"), plot = Norig.ident, device = "pdf", width = 10, height = 10)

  NtagID <- DimPlot(object = obj, group.by = "Sample_Tag", label = T)
  ggsave(filename = paste0(plot_dir, "/", sv.name, "_Sample_Tag.pdf"), plot = NtagID, device = "pdf", width = 10, height = 10)

  NsampelName <- DimPlot(object = obj, group.by = "Sample_Name", label = T)
  ggsave(filename = paste0(plot_dir, "/", sv.name, "_Sample_Name.pdf"), plot = NsampelName, device = "pdf", width = 10, height = 10)

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
    }

    clustre <- clustree(x=obj, prefix="SCT_snn_res.")
    ggsave(filename = paste0(plot_dir, "/", sv.name, "_Clustree.pdf"), plot = clustre, device = "pdf", width = 10, height = 10)
  }

}
