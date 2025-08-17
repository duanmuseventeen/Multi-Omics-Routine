# https://github.com/dpeerlab/Palantir/blob/master/notebooks/Palantir_sample_notebook.ipynb

pip install palantir

import palantir
import scanpy as sc
import pandas as pd
import os

# Plotting
import matplotlib
import matplotlib.pyplot as plt

# warnings
import warnings
from numba.core.errors import NumbaDeprecationWarning

warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
  action="ignore", module="scanpy", message="No data for colormapping"
)

# Inline plotting
%matplotlib inline

ad = sc.read(file_path, backup_url=download_url)
ad

... 

# Palantir next determines the diffusion maps of the data as an estimate of the low dimensional phenotypic manifold of the data.

# Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(ad, n_components=5)

# The low dimensional embeddeing of the data is estimated based on the eigen gap using the following function
ms_data = palantir.utils.determine_multiscale_space(ad)
# If you are specifying the number of eigen vectors manually in the above step, please ensure that the specified parameter is > 2

# I introduce umap from seurat, so don't generate umap again here
# sc.pp.neighbors(ad)
# sc.tl.umap(ad)

# Use scanpy functions to visualize umaps or FDL
sc.pl.embedding(
  ad,
  basis="umap",
  frameon=False,
)

# MAGIC imputation
# MAGIC is an imputation technique developed in the Pe'er lab for single cell data imputation. Palantir uses MAGIC to impute the data for visualization and determining gene expression trends.
imputed_X = palantir.utils.run_magic_imputation(ad)

# sc.pl.embedding(
#   ad,
#   basis="umap",
#   layer="MAGIC_imputed_data",
#   color=["CD34", "MPO", "GATA1", "IRF8"],
#   frameon=False,
# )
# plt.show()

palantir.plot.plot_diffusion_components(ad)
plt.savefig('figures/plot_diffusion_components.pdf')
plt.show()

# Running Palantir
# Palantir can be run by specifying an approxiate early cell.
# 
# Palantir can automatically determine the terminal states as well. In this dataset, we know the terminal states and we will set them using the terminal_states parameter
# 
# The start cell for this dataset was chosen based on high expression of CD34.

# I don't know, so I don't set this para
# terminal_states = pd.Series(
#   ["DC", "Mono", "Ery"],
#   index=["Run5_131097901611291", "Run5_134936662236454", "Run4_200562869397916"],
# )
# The cells can be highlighted on the UMAP map using the highlight_cells_on_umap function

palantir.plot.highlight_cells_on_umap(ad, terminal_states)
plt.show()

# The cell in Neu_CCL4 with the highest CytoTRACE2_Score is selected
start_cell = "sc-07_CCCTCAACATCTGTTT-1"
pr_res = palantir.core.run_palantir(
  ad, start_cell #, num_waypoints=500
)

# Palantir generates the following results
# 
# Pseudotime: Pseudo time ordering of each cell
# Terminal state probabilities: Matrix of cells X terminal states. Each entry represents the probability of the corresponding cell reaching the respective terminal state
# Entropy: A quantiative measure of the differentiation potential of each cell computed as the entropy of the multinomial terminal state probabilities

# Visualizing Palantir results
Palantir results can be visualized on the tSNE or UMAP using the plot_palantir_results function

palantir.plot.plot_palantir_results(ad, s=3)
plt.savefig('figures/plot_palantir_results.pdf')
plt.show()

pr_res = palantir.core.run_palantir(
  ad, # start_cell , num_waypoints=500
)
palantir.plot.plot_palantir_results(ad, s=3)
plt.savefig('figures/plot_palantir_results(nostart).pdf')
plt.show()
# TypeError: run_palantir() missing 1 required positional argument: 'early_cell'

cells = [
  "sc-01_AAACCCAAGCGGGTAT-1",
  "sc-01_AAACCCATCGCAATTG-1",
  "sc-01_AACGAAACAGACAAAT-1",
  "sc-01_AGATGAAGTACCAGAG-1",
]
palantir.plot.plot_terminal_state_probs(ad, cells)
plt.savefig('figures/palantir.plot.plot_terminal_state_probs.pdf')
plt.show()

palantir.plot.highlight_cells_on_umap(ad, ["sc-07_CCCTCAACATCTGTTT-1"])
plt.savefig('figures/early_cell.pdf')
plt.show()
# Gene expression trends
# Gene expression trends over pseudotime provide insights into the dynamic behavior of genes during cellular development or progression. By examining these trends, we can uncover the timing of gene expression changes and identify pivotal regulators of cellular states. Palantir provides tools for computing these gene expression trends.
# Here, we'll outline the steps to compute gene trends over pseudotime using Palantir.
# Selecting cells of a specific trend
# Before computing the gene expression trends, we first need to select cells associated with a specific branch of the pseudotime trajectory. We accomplish this by using the select_branch_cells function. The parameters q and eps are used to control the selection's tolerance. Select small values >=0 to be more sringent and larger values <1 to select more cells.

masks = palantir.presults.select_branch_cells(ad, q=.01, eps=.01)

# Visualizing the branch selection
# Once the cells are selected, it's often helpful to visualize the selection on the pseudotime trajectory to ensure we've isolated the correct cells for our specific trend. We can do this using the plot_branch_selection function:
palantir.plot.plot_branch_selection(ad)
plt.savefig('figures/palantir.plot.plot_branch_selection.pdf')
plt.show()

# To visualize a trajectory on the UMAP, we interpolate the UMAP coordinates of cells specific to each branch across pseudotime, enabling us to draw a continuous path.
palantir.plot.plot_trajectory(ad, "sc-07_AGGCTGCAGCGATGAC-1")
plt.savefig('figures/branch1.pdf')
plt.show()

palantir.plot.plot_trajectory(ad, "sc-04_CGGGTCATCATGGATC-1")
plt.savefig('figures/branch2.pdf')
plt.show()

palantir.plot.plot_trajectory(ad, "sc-01_GAGGCCTGTAGACACG-1")
plt.savefig('figures/branch3.pdf')
plt.show()

# The appearence of this plot is can be adjusted:
  
palantir.plot.plot_trajectory(
  ad, # your anndata
  "DC", # the branch to plot
  cell_color="palantir_entropy", # the ad.obs colum to color the cells by
  n_arrows=10, # the number of arrow heads along the path
  color="red", # the color of the path and arrow heads
  scanpy_kwargs=dict(cmap="viridis"), # arguments passed to scanpy.pl.embedding
  arrowprops=dict(arrowstyle="->,head_length=.5,head_width=.5", lw=3), # appearance of the arrow heads
  lw=3, # thickness of the path
  pseudotime_interval=(0, .9), # interval of the pseudotime to cover with the path
)

# To plot multiple branches use the plot_trajectories with similar customizability. It's default coloring is "palantir_pseudotime" since the overlapping branch masks are hard to vizualize.

palantir.plot.plot_trajectories(ad, pseudotime_interval=(0, .9))
# When using cell_color="branch_selection" be aware of the overlap between branches:
palantir.plot.plot_trajectories(ad, cell_color = "branch_selection", pseudotime_interval=(0, .9))
plt.show()

# Palantir uses Mellon Function Estimator to determine the gene expression trends along different lineages. The marker trends can be determined using the following snippet. This computes the trends for all lineages. A subset of lineages can be used using the lineages parameter.
gene_trends = palantir.presults.compute_gene_trends(
  ad,
  expression_key="MAGIC_imputed_data",
)

genes = ["CD34", "MPO", "GATA1", "IRF8"]
palantir.plot.plot_gene_trends(ad, genes)

palantir.plot.plot_gene_trend_heatmaps(ad, genes)
plt.show()

palantir.plot.plot_trend(ad, "sc-04_CGGGTCATCATGGATC-1", "KLF1", color="nCount_RNA", position_layer="MAGIC_imputed_data")
plt.show()

# Clustering
# Gene expression trends can be clustered and visualized using the following snippet. As an example, the first 1000 genes along the erythroid genes are clustered
mask = (ad.X.toarray() > 1).any(axis=0)
# ad_filtered = ad[:, mask].copy()
more_genes = ad.var_names[mask]
communities = palantir.presults.cluster_gene_trends(ad, "sc-04_CGGGTCATCATGGATC-1", more_genes)
# To achieve the future defaults please pass: flavor="igraph" and n_iterations=2.  directed must also be False to work with igraph's implementation.
# sc.tl.leiden(gt_ad, **kwargs)

palantir.plot.plot_gene_trend_clusters(ad, "sc-04_CGGGTCATCATGGATC-1")
plt.savefig('figures/clustering_branch2.pdf')
plt.show()

# the colname of ad.varm branch are number between 0 and 1, when any two numbers are
# equal, raising errors
del ad.varm['gene_trends_sc-07_AGGCTGCAGCGATGAC-1']

file_path = os.path.join('figures/', "Neu_palantir_processed.h5ad")
ad.write(file_path)