# https://cloud.tencent.com/developer/article/2400875
# https://zhuanlan.zhihu.com/p/647940923

# 1: SeuratDisk ================================================================
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

# 2: scCustomize ===============================================================
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

# 3: sceasy ====================================================================
# https://github.com/meta-cancer/scPLC/blob/main/codes/1b-1.scanpy_seurat2anndata.R

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

# 4: MuDataSeurat ==============================================================
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

# 5: Manual ====================================================================
# https://zhuanlan.zhihu.com/p/25262647576
# 1. 在 Seurat 中导出必要的数据
# 1.1 导出 RNA 计数矩阵
library(Seurat)
library(Matrix)

# 导出 RNA 计数矩阵
counts <- GetAssayData(seurat_object, slot = "counts", assay = "RNA")
writeMM(counts, file = "counts.mtx")

# 导出基因名
write.csv(rownames(counts), file = "genes.csv", row.names = FALSE)

# 导出细胞名
write.csv(colnames(counts), file = "barcodes.csv", row.names = FALSE)
(2) 导出细胞元数据
metadata <- seurat_object@meta.data
write.csv(metadata, "metadata.csv", row.names = TRUE)
(3) 导出降维坐标（如 UMAP, PCA, tSNE, Harmony）
# 这里可以导出不同的降维坐标，如 UMAP, PCA, tSNE, Harmony
umap_coords <- Embeddings(seurat_object, reduction = "umap")
pca_coords <- Embeddings(seurat_object, reduction = "pca")
tsne_coords <- Embeddings(seurat_object, reduction = "tsne")
harmony_coords <- Embeddings(seurat_object, reduction = "harmony")

# 分别保存不同降维坐标
write.csv(umap_coords, "umap_coords.csv", row.names = TRUE)
write.csv(pca_coords, "pca_coords.csv", row.names = TRUE)
write.csv(tsne_coords, "tsne_coords.csv", row.names = TRUE)
write.csv(harmony_coords, "harmony_coords.csv", row.names = TRUE)
2. 在 Scanpy 中转换为 AnnData 格式
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.io import mmread

# 读取计数矩阵
counts = mmread("counts.mtx").T.tocsr()

# 读取基因和细胞名称
genes = pd.read_csv("genes.csv", header=None).squeeze("columns")
cells = pd.read_csv("barcodes.csv", header=None).squeeze("columns")

# **遇到的问题 & 解决思路**
# ❌ 可能出现BUG报错: genes.csv 可能包含无效的第一行，如 'x'
# ✅ 解决：检查第一行并删除
if genes.iloc[0] == "x":
  genes = genes.iloc[1:]
genes = genes[:counts.shape[1]]

# ❌ 可能出现BUG报错: 细胞数不匹配，导致 AnnData 赋值失败
# ✅ 解决：检查长度是否大于计数矩阵，并删除无效行
if len(cells) > counts.shape[0]:
  cells = cells.iloc[1:]
cells = cells[:counts.shape[0]]

# 创建 AnnData 对象
adata = sc.AnnData(X=counts)
adata.var_names = genes
adata.obs_names = cells

# 读取元数据，并确保索引匹配
metadata = pd.read_csv("metadata.csv", index_col=0)
metadata = metadata.loc[adata.obs_names.intersection(metadata.index)]
adata.obs = metadata

# 读取并匹配降维坐标
umap_coords = pd.read_csv("umap_coords.csv", index_col=0)
umap_coords = umap_coords.loc[adata.obs_names]
adata.obsm["X_umap"] = umap_coords.values

harmony_coords = pd.read_csv("harmony_coords.csv", index_col=0)
harmony_coords = harmony_coords.loc[adata.obs_names]
adata.obsm["X_harmony"] = harmony_coords.values

print("AnnData successfully created with Seurat data!")
这样，Seurat 数据成功转换到 Scanpy，同时支持多种降维坐标，并在转换过程中避免了 genes.csv 和 cells.csv 可能引起的错误。

# 3. 统一 Seurat 和 Scanpy 的颜色方案
# 3.1 在 Seurat 中定义颜色并导出
celltype_colors <- c(
  "B_cell" = "#1f78b4",
  "T_cell" = "#33a02c",
  "Monocyte" = "#e31a1c",
  "NK_cell" = "#ff7f00"
)

DimPlot(seurat_object, reduction = "umap", group.by = "cell_type", cols = celltype_colors)

color_mapping <- data.frame(cell_type = names(celltype_colors), color = unname(celltype_colors))
write.csv(color_mapping, "celltype_colors.csv", row.names = FALSE)
# 3.2 在 Scanpy 中应用相同颜色
color_mapping = pd.read_csv("celltype_colors.csv")
color_dict = dict(zip(color_mapping["cell_type"], color_mapping["color"]))
sc.pl.umap(adata, color="cell_type", palette=color_dict)
# 3.3 在 Seurat 中为 orig.ident 定义颜色并导出
origident_colors <- c(
  "Sample1" = "#8c564b",
  "Sample2" = "#e377c2",
  "Sample3" = "#7f7f7f",
  "Sample4" = "#bcbd22"
)
DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident", cols = origident_colors)

origident_mapping <- data.frame(orig.ident = names(origident_colors), color = unname(origident_colors))
write.csv(origident_mapping, "origident_colors.csv", row.names = FALSE)
# 3.4 在 Scanpy 中应用 orig.ident 颜色
origident_mapping = pd.read_csv("origident_colors.csv")
origident_dict = dict(zip(origident_mapping["orig.ident"], origident_mapping["color"]))
sc.pl.umap(adata, color="orig.ident", palette=origident_dict)
