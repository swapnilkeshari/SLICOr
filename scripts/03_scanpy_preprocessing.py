import numpy as np, pandas as pd, scanpy as sc, anndata as ad, matplotlib.pyplot as plt, os, logging, seaborn as sns, pickle
from scipy.stats import median_abs_deviation, hypergeom
import multiprocessing as mp, celloracle as co


wd = '/ocean/projects/cis240075p/skeshari/igvf/bcell2/primaryBCell'
out_path = os.path.join(wd, 'out_data')


oracle = co.load_hdf5(f'{out_path}/intermediate_data/primaryBdata.celloracle.oracle')
links = oracle.get_links(cluster_name_for_GRN_unit="leiden", alpha=10, verbose_level=10)
# Calculate GRN for each population in clustering unit.
links.filter_links() # Getting back to default thresholds (p<0.001 and top 10k) for simulation of GRN -- so that we have enough edges to simulate the GRN
oracle.get_cluster_specific_TFdict_from_Links(links_object=links) #### Runs another round on links.filtered_links to get final links
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
# Saving for KO simulation.
oracle.to_hdf5(f'{out_path}/intermediate_data/oracle_fitted.celloracle.oracle')
links.to_hdf5(f'{out_path}/intermediate_data/oracle_fitted.celloracle.links')
# Saving the GRN for each cluster
links_after_fit = co.Links(name="links_after_fit")
links_after_fit.links_dict = links.links_dict #### Saving (GRN 1 simulation) completely
links_after_fit.filtered_links = links.filtered_links #### Just for initialization
for cluster in links.filtered_links.keys():
    cluster_specific_links = oracle.coef_matrix_per_cluster[cluster].stack().reset_index()
    cluster_specific_links.columns = ['source', 'target', 'coef_mean']
    cluster_specific_links = cluster_specific_links[cluster_specific_links ['coef_mean'] != 0].reset_index(drop=True)
    cluster_specific_links['coef_abs'] = np.abs(cluster_specific_links['coef_mean'])
    cluster_specific_links.to_csv(f'{out_path}/intermediate_data/cluster_{cluster}.csv', index=False)
    links_after_fit.filtered_links[cluster] = cluster_specific_links
# Generating network scores
links_after_fit.get_network_score()
links_after_fit.merged_score.to_csv(f'{out_path}/intermediate_data/network_analysis_scores.csv')

