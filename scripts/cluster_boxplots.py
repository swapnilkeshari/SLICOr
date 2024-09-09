import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import celloracle as co 

######################################################
 # HELPER FUNCTIONS
######################################################
def _get_ix_for_a_cluster(oracle, cluster_column_name, cluster):
    ix = np.arange(oracle.adata.shape[0])[oracle.adata.obs[cluster_column_name] == cluster]
    return ix

######################################################

# Normalize function as provided
def normalize_gradient(gradient, method="sqrt"):
    if method == "sqrt":
        size = np.sqrt(np.power(gradient, 2).sum(axis=1))
        size_sq = np.sqrt(size)
        size_sq[size_sq == 0] = 1
        factor = np.repeat(np.expand_dims(size_sq, axis=1), gradient.shape[1], axis=1)
        normalized_gradient = gradient / factor
    return normalized_gradient

######################################################

def normalize_and_label(diff, oracle, cluster_name):
    """
      Normalize the diffs for each cluster and create labels for plotting.
    """
    cluster_column_name = oracle.cluster_column_name
    # Normalize the diffs for each cluster
    normalized_diffs = [normalize_gradient(diff[_get_ix_for_a_cluster(oracle, cluster_column_name, cluster)])
                        for cluster in cluster_name]

    # Concatenate all normalized diffs into a single array for plotting
    all_normalized_diffs = np.concatenate(normalized_diffs)

    # Since each 'diff' is now a 2D vector, take a norm (e.g., Euclidean norm) to reduce to 1D
    all_normalized_diffs_1d = np.linalg.norm(all_normalized_diffs, axis=1)

    # Create a corresponding list of cluster names for labeling
    cluster_labels = [cluster for cluster in cluster_name for _ in range(len(normalized_diffs[cluster_name.index(cluster)]))]

    return all_normalized_diffs_1d, cluster_labels

######################################################
def _adata_to_color_dict(adata, cluster_use):
    """
    Extract color information from adata and returns as dictionary.

    Args:
        adata (anndata): anndata

        cluster_use (str): column name in anndata.obs

    Returns:
        dictionary: python dictionary, key is cluster name, value is clor name
    """
    color_dict = {}
    for i,j in enumerate(adata.obs[cluster_use].cat.categories):
        color_dict[j] = adata.uns[f"{cluster_use}_colors"][i]
    return color_dict
test = _adata_to_color_dict(oracle.adata, oracle.cluster_column_name)

######################################################
 # MAIN 
######################################################

file_path = 'figures/PRDM1_KO_simulated.celloracle.oracle'

oracle = co.load_hdf5(file_path)

cluster_name = list(oracle.adata.obs["sub_cell_type"])

diff = oracle.delta_embedding - oracle.delta_embedding_random
diff_vec = normalize_and_label(diff = diff, oracle = oracle, cluster_name = cluster_name)
# Now we can plot using seaborn
color_dict = _adata_to_color_dict(oracle.adata, oracle.cluster_column_name)

sns.boxplot(x=diff_vec[1], y= diff_vec[0],palette = color_dict)
plt.title('Normalized Differences in Embedding for each Cluster')
plt.xlabel('Cluster')
plt.ylabel('Normalized Difference in Embedding')
plt.xticks(rotation=90)  # Rotate the x-axis labels for better readability
plt.tight_layout()  # Adjust layout so everything fits without overlapping
plt.show()