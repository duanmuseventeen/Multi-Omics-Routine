#===============================================================================
# Cell Cell Communication with cpdb ============================================
#===============================================================================
# https://www.cellphonedb.org/
# https://github.com/ventolab/CellphoneDB
# https://cloud.tencent.com/developer/article/2518593
# https://github.com/ventolab/CellphoneDB/blob/master/notebooks/T0_DownloadDB.ipynb

# 按照官方提示安装软件即可
# 如果网络不佳，无法下载database文件，可以换用win端下载，完成后复制到linux
# 结果解读请参考官方教程

# cpdb_file_path: (mandatory) path to the database cellphonedb.zip.
# meta_file_path: (mandatory) path to the meta file linking cell barcodes to cluster labels metadata.tsv.
# counts_file_path: (mandatory) paths to normalized counts file (not z-transformed), 
#   either in text format or h5ad (recommended) normalised_log_counts.h5ad.
# microenvs_file_path (optional) path to microenvironment file that groups cell clusters by microenvironments. When providing a microenvironment file, CellphoneDB will restrict the interactions to those cells within a microenvironment.
# active_tf_path: (optional) to the active transcription factors.
# Method 1 ---------------------------------------------------------------------
# R 
Clist <- list()
n <- 1
for (i in unique(seurat.obj$cell_type)) {
  Clist[[n]] <- FindMarkers(seurat.obj, group.by = "cell_type", ident.1 = as.character(i))
  write.csv(Clist[[n]], paste0("cpdb/",i,".csv"))
  
  n <- n + 1
}

normalised_log_counts <- seurat.obj@assays$RNA$data
normalised_log_counts <- as.matrix(normalised_log_counts)
write.csv(normalised_log_counts, "cpdb/normalised_log_counts.csv")
meta.data <- seurat.obj@meta.data
meta.data$barcode <- rownames(meta.data)
write.table(meta.data, "cpdb/bulk_metadata.tsv", sep = "\t")

# python
import pandas as pd
import sys
import os
import anndata

pd.set_option('display.max_columns', 100)
# os.chdir('/home/jovyan/cpdb_tutorial')

cpdb_file_path = '/data/FYM/cpdb/cellphonedb_v500_NatProtocol/v5.0.0/cellphonedb.zip'
meta_file_path = 'bulk_metadata.tsv'
counts_file_path = 'normalised_log_counts.csv'
# microenvs_file_path = 'data/microenvironment.tsv'
out_path = 'method1/'
degs_file_path = 'DEGs.tsv'
# active_tf_path = 'data/active_TFs.tsv'


metadata = pd.read_csv(meta_file_path, sep = '\t')
metadata.head(3)

# Check barcodes in metadata and counts are the same.
# list(adata.obs.index).sort() == list(metadata['barcode_sample']).sort()

# Run basic analysis
# The output of this method will be saved in output_path and also returned to the predefined variables.
from cellphonedb.src.core.methods import cpdb_analysis_method

cpdb_results = cpdb_analysis_method.call(
  cpdb_file_path = cpdb_file_path,           # mandatory: CellphoneDB database zip file.
  meta_file_path = meta_file_path,           # mandatory: tsv file defining barcodes to cell label.
  counts_file_path = counts_file_path,       # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
  counts_data = 'hgnc_symbol',               # defines the gene annotation in counts matrix.
  # microenvs_file_path = microenvs_file_path, # optional (default: None): defines cells per microenvironment.
  score_interactions = True,                 # optional: whether to score interactions or not. 
  output_path = out_path,                    # Path to save results    microenvs_file_path = None,
  separator = '|',                           # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
  threads = 12,                              # number of threads to use in the analysis.
  threshold = 0.2,                           # defines the min % of cells expressing a gene for this to be employed in the analysis.
  result_precision = 3,                      # Sets the rounding for the mean values in significan_means.
  debug = False,                             # Saves all intermediate tables emplyed during the analysis in pkl format.
  output_suffix = None                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
)
# Method 2 (statistical analysis)-----------------------------------------------
Run statistical analysis
The output of this method will be saved in output_path and also returned to the predefined variables.

The statisical method allows the user to downsample the data with the aim of speeding up the results (subsampling arguments). To this end, CellphoneDB employs a geometric sketching procedure (Hie et al. 2019) to preserve the structure of the data without losing information from lowly represented cells. For this tutorial, we have opted to manually downsample the count matrix and the metadata file accordingly.

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
  cpdb_file_path = cpdb_file_path,                 # mandatory: CellphoneDB database zip file.
  meta_file_path = meta_file_path,                 # mandatory: tsv file defining barcodes to cell label.
  counts_file_path = counts_file_path,             # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
  counts_data = 'hgnc_symbol',                     # defines the gene annotation in counts matrix.
  active_tfs_file_path = active_tf_path,           # optional: defines cell types and their active TFs.
  microenvs_file_path = microenvs_file_path,       # optional (default: None): defines cells per microenvironment.
  score_interactions = True,                       # optional: whether to score interactions or not. 
  iterations = 1000,                               # denotes the number of shufflings performed in the analysis.
  threshold = 0.2,                                 # defines the min % of cells expressing a gene for this to be employed in the analysis.
  threads = 5,                                     # number of threads to use in the analysis.
  debug_seed = 42,                                 # debug randome seed. To disable >=0.
  result_precision = 3,                            # Sets the rounding for the mean values in significan_means.
  pvalue = 0.05,                                   # P-value threshold to employ for significance.
  subsampling = False,                             # To enable subsampling the data (geometri sketching).
  subsampling_log = False,                         # (mandatory) enable subsampling log1p for non log-transformed data inputs.
  subsampling_num_pc = 100,                        # Number of componets to subsample via geometric skectching (dafault: 100).
  subsampling_num_cells = 1000,                    # Number of cells to subsample (integer) (default: 1/3 of the dataset).
  separator = '|',                                 # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
  debug = False,                                   # Saves all intermediate tables employed during the analysis in pkl format.
  output_path = out_path,                          # Path to save results.
  output_suffix = None                             # Replaces the timestamp in the output files by a user defined string in the  (default: None).
)

# Method 3 (differential expression)--------------------------------------------
Run CellphoneDB with differential analysis (method 3)
The output of this method will be saved in output_path and also assigned to the predefined variables.

from cellphonedb.src.core.methods import cpdb_degs_analysis_method

cpdb_results = cpdb_degs_analysis_method.call(
  cpdb_file_path = cpdb_file_path,                            # mandatory: CellphoneDB database zip file.
  meta_file_path = meta_file_path,                            # mandatory: tsv file defining barcodes to cell label.
  counts_file_path = counts_file_path,                        # mandatory: normalized count matrix - a path to the counts file, or an in-memory AnnData object
  degs_file_path = degs_file_path,                            # mandatory: tsv file with DEG to account.
  counts_data = 'hgnc_symbol',                                # defines the gene annotation in counts matrix.
  active_tfs_file_path = active_tf_path,                      # optional: defines cell types and their active TFs.
  microenvs_file_path = microenvs_file_path,                  # optional (default: None): defines cells per microenvironment.
  score_interactions = True,                                  # optional: whether to score interactions or not. 
  threshold = 0.2,                                            # defines the min % of cells expressing a gene for this to be employed in the analysis.
  result_precision = 3,                                       # Sets the rounding for the mean values in significan_means.
  separator = '|',                                            # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
  debug = False,                                              # Saves all intermediate tables emplyed during the analysis in pkl format.
  output_path = out_path,                                     # Path to save results
  output_suffix = None,                                       # Replaces the timestamp in the output files by a user defined string in the  (default: None)
  threads = 25
)







