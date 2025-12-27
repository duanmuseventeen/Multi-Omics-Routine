# Reference---------------------------------------------------------------------
# https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html

# Getting started with anndata
# Authors: Adam Gayoso, Alex Wolf

# Note
# This tutorial is based on a blog posts by Adam in 2021 and Alex in 2017.

# In this tutorial, we introduce basic properties of the central object, AnnData (“Annotated Data”).
# 
# AnnData is specifically designed for matrix-like data. By this we mean that we have 
# observations, each of which can be represented as 
# -dimensional vectors, where each dimension corresponds to a variable or feature. Both the rows and columns of this 
# matrix are special in the sense that they are indexed.

# For instance, in scRNA-seq data, each row corresponds to a cell with a barcode, and each column corresponds to a gene with a gene id. Furthermore, for each cell and each gene we might have additional metadata, like (1) donor information for each cell, or (2) alternative gene symbols for each gene. Finally, we might have other unstructured metadata like color palletes to use for plotting. Without going into every fancy Python-based data structure, we think that still today no other alternative really exists that:
#   
# · Handles sparsity
# · Handles unstructured data
# · Handles observation- and feature-level metadata
# · Is user-friendly

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

# Initializing AnnData----------------------------------------------------------
# Let’s start by building a basic AnnData object with some sparse count information, perhaps representing gene expression counts.

counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)
adata
AnnData object with n_obs × n_vars = 100 × 2000
We see that AnnData provides a representation with summary stastics of the data The initial data we passed are accessible as a sparse matrix using adata.X.

adata.X
<100x2000 sparse matrix of type '<class 'numpy.float32'>'
with 126526 stored elements in Compressed Sparse Row format>
  Now, we provide the index to both the obs and var axes using .obs_names (resp. .var_names).

adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]
print(adata.obs_names[:10])
Index(['Cell_0', 'Cell_1', 'Cell_2', 'Cell_3', 'Cell_4', 'Cell_5', 'Cell_6',
       'Cell_7', 'Cell_8', 'Cell_9'],
      dtype='object')
# Subsetting AnnData------------------------------------------------------------
# These index values can be used to subset the AnnData, which provides a view of the AnnData object. We can imagine this to be useful to subset the AnnData to particular cell types or gene modules of interest. The rules for subsetting AnnData are quite similar to that of a Pandas DataFrame. You can use values in the obs/var_names, boolean masks, or cell index integers.

adata[["Cell_1", "Cell_10"], ["Gene_5", "Gene_1900"]]
View of AnnData object with n_obs × n_vars = 2 × 2
# Adding aligned metadata-------------------------------------------------------
# Observation/Variable level
# So we have the core of our object and now we’d like to add metadata at both the observation and variable levels. This is pretty simple with AnnData, both adata.obs and adata.var are Pandas DataFrames.

ct = np.random.choice(["B", "T", "Monocyte"], size=(adata.n_obs,))
adata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency
adata.obs
cell_type
Cell_0	B
Cell_1	B
Cell_2	B
Cell_3	Monocyte
Cell_4	Monocyte
...	...
Cell_95	Monocyte
Cell_96	B
Cell_97	Monocyte
Cell_98	B
Cell_99	T
100 rows × 1 columns

We can also see now that the AnnData representation has been updated:
  
adata
# AnnData object with n_obs × n_vars = 100 × 2000
# obs: 'cell_type'
# Subsetting using metadata-----------------------------------------------------
We can also subset the AnnData using these randomly generated cell types:
  
bdata = adata[adata.obs.cell_type == "B"]
bdata
# View of AnnData object with n_obs × n_vars = 26 × 2000
# obs: 'cell_type'
# Observation/variable-level matrices-------------------------------------------
We might also have metadata at either level that has many dimensions to it, such as a UMAP embedding of the data. For this type of metadata, AnnData has the .obsm/.varm attributes. We use keys to identify the different matrices we insert. The restriction of .obsm/.varm are that .obsm matrices must length equal to the number of observations as .n_obs and .varm matrices must length equal to .n_vars. They can each independently have different number of dimensions.

Let’s start with a randomly generated matrix that we can interpret as a UMAP embedding of the data we’d like to store, as well as some random gene-level metadata:
  
  adata.obsm["X_umap"] = np.random.normal(0, 1, size=(adata.n_obs, 2))
adata.varm["gene_stuff"] = np.random.normal(0, 1, size=(adata.n_vars, 5))
adata.obsm
AxisArrays with keys: X_umap
Again, the AnnData representation is updated.

adata
AnnData object with n_obs × n_vars = 100 × 2000
obs: 'cell_type'
obsm: 'X_umap'
varm: 'gene_stuff'
A few more notes about .obsm/.varm

The “array-like” metadata can originate from a Pandas DataFrame, scipy sparse matrix, or numpy dense array.

When using scanpy, their values (columns) are not easily plotted, where instead items from .obs are easily plotted on, e.g., UMAP plots.

# Unstructured metadata---------------------------------------------------------
AnnData has .uns, which allows for any unstructured metadata. This can be anything, like a list or a dictionary with some general information that was useful in the analysis of our data.

adata.uns["random"] = [1, 2, 3]
adata.uns
OverloadedDict, wrapping:
  OrderedDict([('random', [1, 2, 3])])
With overloaded keys:
  ['neighbors'].
# Layers------------------------------------------------------------------------
# Finally, we may have different forms of our original core data, perhaps one that is normalized and one that is not. These can be stored in different layers in AnnData. For example, let’s log transform the original data and store it in a layer:
  
adata.layers["log_transformed"] = np.log1p(adata.X)
adata
# AnnData object with n_obs × n_vars = 100 × 2000
# obs: 'cell_type'
# uns: 'random'
# obsm: 'X_umap'
# varm: 'gene_stuff'
# layers: 'log_transformed'
# Conversion to DataFrames------------------------------------------------------
We can also ask AnnData to return us a DataFrame from one of the layers:
  
adata.to_df(layer="log_transformed")
# Gene_0	Gene_1	Gene_2	Gene_3	Gene_4	Gene_5	Gene_6	Gene_7	Gene_8	Gene_9	...	Gene_1990	Gene_1991	Gene_1992	Gene_1993	Gene_1994	Gene_1995	Gene_1996	Gene_1997	Gene_1998	Gene_1999
# Cell_0	1.098612	0.693147	0.000000	0.693147	0.693147	0.000000	1.386294	0.693147	0.693147	0.693147	...	1.098612	0.000000	0.693147	0.000000	0.000000	0.693147	0.693147	0.000000	1.098612	0.693147
# Cell_1	0.000000	1.098612	0.693147	0.000000	0.693147	0.693147	0.693147	0.693147	0.693147	0.000000	...	0.693147	0.000000	0.000000	0.000000	0.693147	1.098612	1.098612	0.000000	0.000000	1.386294
# Cell_2	0.693147	0.693147	0.000000	0.693147	1.098612	0.693147	0.693147	0.000000	0.693147	1.098612	...	0.000000	0.000000	0.693147	0.693147	1.386294	0.693147	1.098612	0.000000	0.000000	0.000000
# Cell_3	0.000000	1.098612	0.000000	0.693147	1.791759	0.693147	0.000000	0.000000	1.098612	0.000000	...	1.609438	1.098612	0.693147	0.000000	1.098612	0.000000	0.693147	0.693147	0.693147	0.693147
# Cell_4	0.693147	0.000000	0.693147	0.000000	0.693147	0.693147	0.000000	0.693147	0.693147	1.098612	...	0.693147	1.098612	0.000000	0.000000	0.000000	1.098612	0.000000	1.098612	1.609438	0.693147
# ...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...	...
# Cell_95	0.693147	0.693147	0.000000	1.098612	0.693147	1.098612	0.000000	0.000000	0.000000	0.000000	...	0.000000	0.000000	0.000000	1.098612	0.000000	0.000000	0.693147	0.693147	0.693147	0.000000
# Cell_96	0.693147	1.098612	1.386294	0.693147	1.098612	0.000000	1.609438	0.693147	0.693147	0.693147	...	1.098612	0.000000	0.000000	0.000000	0.000000	0.693147	1.386294	0.000000	1.386294	0.000000
# Cell_97	0.000000	0.693147	1.098612	0.693147	0.000000	0.000000	0.000000	0.000000	1.386294	0.000000	...	0.693147	0.000000	0.693147	0.000000	1.386294	1.386294	0.000000	0.000000	0.693147	0.000000
# Cell_98	0.693147	1.098612	0.000000	0.693147	0.693147	0.000000	0.693147	0.000000	0.693147	1.098612	...	0.693147	0.000000	0.000000	1.098612	0.693147	0.000000	0.693147	0.693147	0.693147	1.098612
# Cell_99	0.693147	0.693147	0.000000	1.791759	0.000000	1.098612	0.000000	1.098612	1.609438	0.693147	...	1.098612	1.098612	0.693147	0.693147	1.098612	0.693147	0.000000	0.000000	0.693147	0.693147
# 100 rows × 2000 columns

We see that the .obs_names/.var_names are used in the creation of this Pandas object.

# Writing the results to disk---------------------------------------------------
# AnnData comes with its own persistent HDF5-based file format: h5ad. If string columns with small number of categories aren’t yet categoricals, AnnData will auto-transform to categoricals.

adata.write('my_results.h5ad', compression="gzip")
!h5ls 'my_results.h5ad'
X                        Group
layers                   Group
obs                      Group
obsm                     Group
obsp                     Group
uns                      Group
var                      Group
varm                     Group
varp                     Group
# Wrapping up the introduction--------------------------------------------------
AnnData has become the standard for single-cell analysis in Python and for good reason – it’s straightforward to use and faciliatates more reproducible analyses with it’s key-based storage. It’s even becoming easier to convert to the popular R-based formats for single-cell analysis.

Keep reading on to better understand “views”, on-disk backing, and other details.

# Views and copies
For the fun of it, let’s look at another metadata use case. Imagine that the observations come from instruments characterizing 10 readouts in a multi-year study with samples taken from different subjects at different sites. We’d typically get that information in some format and then store it in a DataFrame:
  
  obs_meta = pd.DataFrame({
    'time_yr': np.random.choice([0, 2, 4, 8], adata.n_obs),
    'subject_id': np.random.choice(['subject 1', 'subject 2', 'subject 4', 'subject 8'], adata.n_obs),
    'instrument_type': np.random.choice(['type a', 'type b'], adata.n_obs),
    'site': np.random.choice(['site x', 'site y'], adata.n_obs),
  },
  index=adata.obs.index,    # these are the same IDs of observations as above!
  )
This is how we join the readout data with the metadata. Of course, the first argument of the following call for X could also just be a DataFrame.

adata = ad.AnnData(adata.X, obs=obs_meta, var=adata.var)
Now we again have a single data container that keeps track of everything.

print(adata)
AnnData object with n_obs × n_vars = 100 × 2000
obs: 'time_yr', 'subject_id', 'instrument_type', 'site'
Subsetting the joint data matrix can be important to focus on subsets of variables or observations, or to define train-test splits for a machine learning model.

# Note
# Similar to numpy arrays, AnnData objects can either hold actual data or reference another AnnData object. In the later case, they are referred to as “view”.

Subsetting AnnData objects always returns views, which has two advantages:
  
  no new memory is allocated

it is possible to modify the underlying AnnData object

You can get an actual AnnData object from a view by calling .copy() on the view. Usually, this is not necessary, as any modification of elements of a view (calling .[] on an attribute of the view) internally calls .copy() and makes the view an AnnData object that holds actual data. See the example below.

adata
# AnnData object with n_obs × n_vars = 100 × 2000
# obs: 'time_yr', 'subject_id', 'instrument_type', 'site'
# Get access to the first 5 rows for two variables.

# Note
# Indexing into AnnData will assume that integer arguments to [] behave like .iloc in pandas, whereas string arguments behave like .loc. AnnData always assumes string indices.

adata[:5, ['Gene_1', 'Gene_3']]
View of AnnData object with n_obs × n_vars = 5 × 2
obs: 'time_yr', 'subject_id', 'instrument_type', 'site'
This is a view! If we want an AnnData that holds the data in memory, let’s call .copy()

adata_subset = adata[:5, ['Gene_1', 'Gene_3']].copy()
For a view, we can also set the first 3 elements of a column.

print(adata[:3, 'Gene_1'].X.toarray().tolist())
adata[:3, 'Gene_1'].X = [0, 0, 0]
print(adata[:3, 'Gene_1'].X.toarray().tolist())
[[1.0], [2.0], [1.0]]
[[0.0], [0.0], [0.0]]
If you try to access parts of a view of an AnnData, the content will be auto-copied and a data-storing object will be generated.

adata_subset = adata[:3, ['Gene_1', 'Gene_2']]
adata_subset
View of AnnData object with n_obs × n_vars = 3 × 2
obs: 'time_yr', 'subject_id', 'instrument_type', 'site'
adata_subset.obs['foo'] = range(3)
# /var/folders/bd/43q20k0n6z15tdfzxvd22r7c0000gn/T/ipykernel_25768/2955902014.py:1: ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
adata_subset.obs['foo'] = range(3)
Now adata_subset stores the actual data and is no longer just a reference to adata.

adata_subset
AnnData object with n_obs × n_vars = 3 × 2
obs: 'time_yr', 'subject_id', 'instrument_type', 'site', 'foo'
Evidently, you can use all of pandas to slice with sequences or boolean indices.

adata[adata.obs.time_yr.isin([2, 4])].obs.head()
time_yr	subject_id	instrument_type	site
Cell_1	4	subject 4	type b	site y
Cell_2	4	subject 1	type a	site y
Cell_3	4	subject 1	type b	site x
Cell_4	2	subject 1	type b	site x
Cell_6	2	subject 1	type b	site x
# Partial reading of large data-------------------------------------------------
# If a single .h5ad is very large, you can partially read it into memory by using backed mode:
  
adata = ad.read('my_results.h5ad', backed='r')
adata.isbacked
# True
# If you do this, you’ll need to remember that the AnnData object has an open connection to the file used for reading:
  
adata.filename
# PosixPath('my_results.h5ad')
As we’re using it in read-only mode, we can’t damage anything. To proceed with this tutorial, we still need to explicitly close it:
  
adata.file.close()
# As usual, you should rather use with statements to avoid dangling open files (up-coming feature).

# Manipulating the object on disk is possible, but experimental for sparse data. Hence, we leave it out of this tutorial.