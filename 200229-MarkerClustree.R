##This function makes a clustree, where the nodes are colored according to the mean expression of a gene of interest.
#@Param seurat: Give a seurat object to use for clustree
#@Param gene: Give gene of interest. If multiple genes are given the function will output a list of all plots
library(clustree)
library(Seurat)
MarkerClustree <- function(seurat, gene){
  #if length of gene is 1 (i.e. only one gene is given), then only ouput one plot. Otherwise make a list of plots
  if(length(gene)==1){
    #Make a data.frame with meta information from seurat object and gene of interest expression, then format column names
    plot <- clustree(x = seurat, prefix = "RNA_snn_res.", node_colour=gene, node_colour_aggr = "median")
    return(plot)

  } else{
    plot_list <- list() #Output list for final plots
    for(gen in gene){
      #Create clustree and save output in list
      plot <- clustree(x = seurat, prefix = "RNA_snn_res.", node_colour=gen, node_colour_aggr = "median")
      plot_list[[gen]] <- plot
    }
    return(plot_list)
  }
}
