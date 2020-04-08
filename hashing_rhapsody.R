##Function to add sample demultiplexing information to a Rhapsody Seurat object.
hashing_rhapsody <- function(rhap_obj, names_given = FALSE){
  #read in hashing file
  hashing_info <- read.csv("/share/ScratchGeneral/jerven/Rhapsody/Nova/BothLanes/_1_AAGAGGCA_Sample_Tag_Calls.csv", comment.char="#", stringsAsFactors=FALSE)

  #make vectors containing tag ID and sample name
  cell_ids <- hashing_info$Cell_index
  cell_tag <- hashing_info$Sample_Tag
  cell_sampleName <- hashing_info$Sample_Name
  #set vector names to the cell IDs
  names(cell_tag) <- cell_ids
  names(cell_sampleName) <- cell_ids
  #Add meta data to rhapsody object
  rhap_obj <- AddMetaData(object = rhap_obj, metadata = cell_tag, col.name = "Sample_Tag")

  #If sample names were not give, translate the tag IDs into the right sample names
  if(!names_given){
    tag_translation <- data.frame("tag"=c("SampleTag01_mm", "SampleTag03_mm", "SampleTag02_mm", "Undetermined", "Multiplet", "SampleTag05_mm"),
                                  "name"=c("4T1_tumor", "67NR_tumor", "4T1_lung", "Undetermined", "Multiplet", "67NR_lung"))
    for(tag in tag_translation$tag){
    cell_sampleName[cell_sampleName == tag] <- as.character(tag_translation[tag_translation$tag == tag, "name"])
    }
  }

  #Add meta data of sample names to rhapsody object
  rhap_obj <- AddMetaData(object = rhap_obj, metadata = cell_sampleName, col.name = "Sample_Name")

  return(rhap_obj)
}
