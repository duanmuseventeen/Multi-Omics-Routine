# https://cloud.tencent.com/developer/article/2400875
# https://zhuanlan.zhihu.com/p/647940923

# 1: SeuratDisk 
# https://mojaveazure.github.io/seurat-disk/articles/convert-anndata.html
require(SeuratDisk)
system.time({
  SaveH5Seurat(seurat.obj, filename = "G:/Seurat2h5/T_cell")
  Convert("G:/Seurat2h5/T_cell.h5seurat", dest = "h5ad")
})
if(FALSE){
  Creating h5Seurat file for version 3.1.5.9900
  Adding cell embeddings for pca
  Adding loadings for pca
  No projected loadings for pca
  Adding standard deviations for pca
  No JackStraw data for pca
  Adding cell embeddings for harmony
  Adding loadings for harmony
  Adding projected loadings for harmony
  Adding standard deviations for harmony
  No JackStraw data for harmony
  Adding cell embeddings for umap
  No loadings for umap
  No projected loadings for umap
  No standard deviations for umap
  No JackStraw data for umap
  Validating h5Seurat file
  Adding data from RNA as X
  错误于assay.group$obj_copy_to(dst_loc = dfile, dst_name = "X", src_name = x.data): 
    HDF5-API Errors:
    error #000: ../../src/H5Ocopy.c in H5Ocopy(): line 240: unable to copy object
  class: HDF5
  major: Object header
  minor: Unable to copy object
  
  error #001: ../../src/H5VLcallback.c in H5VL_object_copy(): line 5495: object copy failed
  class: HDF5
  major: Virtual Object Layer
  minor: Unable to copy object
  
  error #002: ../../src/H5VLcallback.c in H5VL__object_copy(): line 5456: object copy failed
  class: HDF5
  major: Virtual Object Layer
  minor: Unable to copy object
  
  error #003: ../../src/H5VLnative_object.c in H5VL__native_object_copy(): line 125: unable to copy object
  class: HDF5
  major: Object header
  minor: Unable to copy object
  
  error #004: ../../src/H5Ocopy.c in H5O__copy(): line 291: source object not found
  class: HDF5
  major: Symbol table
  minor: Object not found
  
  error #005: ../../src/H5Gloc.c in H5G_loc_find(): line 442: 
  Timing stopped at: 18.71 1.31 23.46
}

# 2: scCustomize 
# https://github.com/satijalab/seurat/discussions/8642
require(scCustomize)
system.time({
  as.anndata(x = seurat.obj, 
             file_path = "G:/Seurat2h5/", 
             file_name = "T_cell(scCustomize).h5ad", 
             assay = "RNA", 
             main_layer = "data", 
             other_layers = c("counts"), 
             transfer_dimreduc = TRUE, 
             verbose = TRUE)
})
if(FALSE){
  • Checking Seurat object validity
  • Extracting Data from RNA assay to transfer to anndata.
  The following columns were removed as they contain identical values for all rows:
    ℹ percent.mt and scDblFinder.class
  • Creating anndata object.
  • Writing anndata file: "G:\\Seurat2h5/T_cell(scCustomize).h5ad"
  用户  系统  流逝 
  18.84  2.02 51.36 
  警告信息:
    Adding a command log without an assay associated with it 
}

# 3: sceasy
# BiocManager::install(c("LoomExperiment"))
# devtools::install_github("cellgeni/sceasy")
library(sceasy)
library(reticulate)
# use_condaenv('EnvironmentName')
loompy <- reticulate::import('loompy')

# https://github.com/cellgeni/sceasy/issues/82
# for seurat v5, run this code to downgrade the assay:
seurat.obj[["RNA"]] <- as(seurat.obj[["RNA"]], "Assay")
system.time({
  sceasy::convertFormat(seurat.obj, 
                        from="seurat", 
                        to="anndata",
                        outFile='G:/Seurat2h5/T_cell(sceasy).h5ad')
})
if(FALSE){
  用户  系统  流逝 
  11.30  1.12 22.33 
  警告信息:
    In .regularise_df(obj@meta.data, drop_single_values = drop_single_values) :
    Dropping single category variables:percent.mt, scDblFinder.class
}

# 4: MuDataSeurat
# omit

# 5: dior
# devtools::install_github('JiekaiLab/dior')
# library(dior)
system.time({
  dior::write_h5(seurat.obj, 
                 file='G:/Seurat2h5/T_cell(dior).h5', 
                 object.type = 'seurat', 
                 assay.name = "RNA")
})
if(FALSE){
  用户  系统  流逝 
  11.98  0.59 20.39 
}