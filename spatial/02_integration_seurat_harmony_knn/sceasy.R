library(sceasy)
library(reticulate)

# use the specified Conda environment
use_condaenv('sceasy')
loompy <- reticulate::import('loompy')

convert_seurat_to_anndata <- function(seurat_object, output_path) {
  # check the input parameters
  if (missing(seurat_object) || missing(output_path)) {
    stop("Please provide the Seurat object and the output path.")
  }
  
  # execute the format conversion
  sceasy::convertFormat(seurat_object, from = "seurat", to = "anndata", outFile = output_path)
  
  message("Conversion successful, the output file is in: ", output_path)
}

# example usage
# convert_seurat_to_anndata(seurat_object, "output_file.h5ad")

convert_anndata_to_seurat <- function(anndata_file, output_path) {
  # check the input parameters
  if (missing(anndata_file) || missing(output_path)) {
    stop("Please provide the Anndata object and the output path.")
  }
  
  # execute the format conversion
  sceasy::convertFormat(anndata_file, from = "anndata", to = "seurat", outFile = output_path)
  
  message("Conversion successful, the output file is in: ", output_path)
}

