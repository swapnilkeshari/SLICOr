import os, math
import sys, itertools
from scipy.stats import hypergeom
import multiprocessing
from multiprocessing.pool import Pool
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import networkx as nx
import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',level=logging.INFO,datefmt='%d-%b-%y %H:%M:%S')
import celloracle as co
logging.info(co.__version__)

def _get_ix_for_a_cluster(oracle, cluster_column_name, cluster):
    ix = np.arange(oracle.adata.shape[0])[oracle.adata.obs[cluster_column_name] == cluster]
    return ix

# Normalize function as provided
def normalize_gradient(gradient, method="sqrt"):
    if method == "sqrt":
        size = np.sqrt(np.power(gradient, 2).sum(axis=1))
        size_sq = np.sqrt(size)
        size_sq[size_sq == 0] = 1
        factor = np.repeat(np.expand_dims(size_sq, axis=1), gradient.shape[1], axis=1)
        normalized_gradient = gradient / factor
    return normalized_gradient

def normalize_and_label(diff, oracle, cluster_name,cluster_column_name):
    """
      Normalize the diffs for each cluster and create labels for plotting.
    """
    cluster_column_name = cluster_column_name
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
def _adata_to_color_dict(adata, cluster_use):
    color_dict = {}
    for i,j in enumerate(adata.obs[cluster_use].cat.categories):
        color_dict[j] = adata.uns[f"{cluster_use}_colors"][i]
    return color_dict

def ko_simulation(oracle,goi,cluster_column_name="sub_cell_type"):
    try:
        print(f'Process {multiprocessing.current_process().name} started working on gene {goi}', flush=True)
        print(multiprocessing.current_process(),flush=True)
        oracle.simulate_shift(perturb_condition={goi: 0.0},n_propagation=3)
        # Get transition probability
        oracle.estimate_transition_prob(n_neighbors=200,knn_random=True, sampled_fraction=1)

        # Calculate embedding 
        oracle.calculate_embedding_shift(sigma_corr=0.05)
        print("Calculated the embedding shift",flush=True)
        cluster_column_name="sub_cell_type"
        cluster_name = list(oracle.adata.obs[cluster_column_name].unique())
        diff = oracle.delta_embedding - oracle.delta_embedding_random
        diff_vec = pd.DataFrame(normalize_and_label(diff = diff, oracle = oracle, cluster_name = cluster_name,cluster_column_name = cluster_column_name)).T
        # Now we can plot using seaborn
        l_interest=['5_Naive','3_GC','1_ABC','7_PB','6_MBC from GC']
        diff_vec=diff_vec[diff_vec[1].isin(l_interest)]
        diff_vec.to_csv(f"/ix/djishnu/Swapnil/CellOracle/primaryBCell/out_data/nick_cshl_apr_2024/{goi}_delta_embedding_KO.csv")

        color_dict = _adata_to_color_dict(oracle.adata, cluster_column_name)
        sns.boxplot(x=diff_vec[1], y= diff_vec[0],palette = color_dict)
        plt.title(f'Normalized Differences in Embedding for each Cluster for {goi}')
        plt.xlabel('Cluster')
        plt.ylabel('Normalized Difference in Embedding')
        plt.xticks(rotation=90)  # Rotate the x-axis labels for better readability
        plt.tight_layout()  # Adjust layout so everything fits without overlapping
        plt.savefig(f"figures/delta_embedding:{goi}_KO.jpg")
        print(f'Process {multiprocessing.current_process().name} ended working on gene {goi}', flush=True)
    except Exception as e:
        print(f'Process {multiprocessing.current_process().name} failed for gene {goi} with error {e}', flush=True)
        return None
    return oracle

if __name__ == "__main__":
    plt.rcParams['figure.figsize'] = [6, 4.5]
    plt.rcParams["savefig.dpi"] = 300

    save_folder = "../out_data/figures"
    os.makedirs(save_folder, exist_ok=True)
    oracle = co.load_hdf5("../out_data/primaryBdata.celloracle.oracle")
    print(oracle)
    links = co.load_hdf5("../out_data/links.celloracle.links")

    links.filter_links()
    oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

    # %%
    base_GRN=pd.read_csv('../out_data/base_GRN_dataframe.csv')
    maybe_CO_TFs=list(base_GRN.groupby(['gene_short_name']).sum().columns)
    list_of_LFs = pd.read_csv("/ix/djishnu/Zarifeh/B_cell_oracle/GC_PB_Multiome/LFs.txt",sep='\t')
    # Significant LF Paths
    sig_paths=list(list_of_LFs['Sig_LFs'])
    er_sig_sa=[]
    for path in sig_paths:
        if type(path) == str:
            er_sig_sa=er_sig_sa+list(pd.read_csv(path,sep='\t')['names'])
    # Interacting LF Paths
    interacting_paths=list(list_of_LFs['Interacting_LFs'])
    er_sig_interacting=[]
    for path in interacting_paths:
        if type(path) == str:
            er_sig_interacting=er_sig_interacting+list(pd.read_csv(path,sep='\t')['names'])

    set_er_sig_sa = set(er_sig_sa)
    set_maybe_CO_TFs = set(maybe_CO_TFs)
    set_er_sig_interacting = set(er_sig_interacting)

    er_TFs = (set_er_sig_sa & set_maybe_CO_TFs) | (set_er_sig_interacting & set_maybe_CO_TFs)
    er_genes = (set_er_sig_sa - set_maybe_CO_TFs) | (set_er_sig_interacting - set_maybe_CO_TFs)
    cluster_list = ['3','7']
    TF_er_genes,not_in_graph=[],[]
    min_weight_threshold = 0.5
    total_er_genes_in_LFs = len(er_TFs)+len(er_genes)
    total_er_genes_in_start = 6700
    # Initialize the dictionary with keys and empty values
    cluster_graphs = {key: None for key in cluster_list}

    for cluster in cluster_list:
        source=list(links.filtered_links[cluster]['source'])
        target=list(links.filtered_links[cluster]['target'])
        weight=list(links.filtered_links[cluster]['-logp']*links.filtered_links[cluster]['coef_mean'])
        edge_data = zip(source, target, weight)
        filtered_edge_data = [(src, tgt, wt) for src, tgt, wt in edge_data if abs(wt) >= min_weight_threshold]
        G = nx.Graph()
        G.add_weighted_edges_from(filtered_edge_data)
        cluster_graphs[cluster] = G
        for gene in er_genes:
            if gene in G.nodes:
                TF_er_genes.extend(set(cluster_graphs[cluster][gene].keys()) & set(maybe_CO_TFs))
    TF_er_genes = list(set(TF_er_genes))

    def enrichment_dict(TF_list,cluster_list,cluster_graphs,er_genes,total_er_genes_in_LFs,total_er_genes_in_start):
        TF_dict = {key: {cluster: (None,None) for cluster in cluster_list} for key in TF_list}
        for TF in TF_list:
            for cluster in cluster_list:
                if TF in cluster_graphs[cluster].nodes:
                    common=len(list(set(cluster_graphs[cluster][TF].keys()).intersection(set(er_genes))))
                    dwngene=len(list(cluster_graphs[cluster][TF].keys())) #TF_sub_net_df['target']
                    if dwngene==0 or common==0:
                        TF_dict[TF][cluster]= (0,1)
                        # logging.warn(f'No downstream {dwngene} or common {common} genes found for {TF} in cluster {cluster}\n')
                    else:
                        M = total_er_genes_in_start  # Total population size
                        n = dwngene   # Total number of successes in the population
                        N = total_er_genes_in_LFs   # Sample size
                        X = common   # Observed number of successes in the sample
                        p_value = 1 - hypergeom.cdf(X-1, M, n, N)   # Compute the p-value for observing X or more successes
                        TF_dict[TF][cluster] = (math.log2((common/(dwngene))/(total_er_genes_in_LFs/total_er_genes_in_start)), p_value)
        return TF_dict

    er_TF_dict = enrichment_dict(er_TFs,cluster_list,cluster_graphs,er_genes,total_er_genes_in_LFs,total_er_genes_in_start)
    TF_er_genes_dict = enrichment_dict(TF_er_genes,cluster_list,cluster_graphs,er_genes,total_er_genes_in_LFs,total_er_genes_in_start)

    pc_er_TFs = list(itertools.combinations(er_TFs, 2))
    pc_TF_er_genes = list(itertools.combinations(TF_er_genes, 2))
    pc_TF_all = list(itertools.combinations(er_TFs.union(TF_er_genes), 2))

    TFs_dict = {key: {cluster: (None,None) for cluster in cluster_list} for key in pc_TF_all}
    for TFs in pc_TF_all:
        for cluster in cluster_list:
            if TFs[0] in cluster_graphs[cluster].nodes and TFs[1] in cluster_graphs[cluster].nodes:
                common_tf1=set(cluster_graphs[cluster][TFs[0]].keys()).intersection(set(er_genes))
                common_tf2=set(cluster_graphs[cluster][TFs[1]].keys()).intersection(set(er_genes))
                common=len(common_tf1.intersection(common_tf2))
                dwngene_1=set(cluster_graphs[cluster][TFs[0]].keys())
                dwngene_2=set(cluster_graphs[cluster][TFs[1]].keys())
                dwngene=len(list(dwngene_1.intersection(dwngene_2)))
                if dwngene==0 or common==0:
                    TFs_dict[TFs][cluster]=(0,1)
                    # logging.warn(f'No downstream {dwngene} or common {common} genes found for {TF} in cluster {cluster}\n')
                else:
                    M = total_er_genes_in_start  # Total population size
                    n = dwngene   # Total number of successes in the population
                    N = total_er_genes_in_LFs   # Sample size
                    X = common   # Observed number of successes in the sample
                    p_value = 1 - hypergeom.cdf(X-1, M, n, N)   # Compute the p-value for observing X or more successes
                    TFs_dict[TFs][cluster]=(math.log2((common/(dwngene))/(total_er_genes_in_LFs/total_er_genes_in_start)),p_value)

    # %%
    with Pool(processes=8) as pool:
        tf_of_interest = ["BATF","BATF3","BCL6","EGR1","FOS","IKZF1","IRF1","IRF4","IRF8","JUNB","JUND","MEF2A","MEF2C","MYB","NFATC1","NFATC2","NFIL3","NFKB1","NFKB2","PAX5","PRDM1","RUNX1","SP3","SPI1","SPIB","STAT1","STAT5A","TCF12","VDR","XBP1","ZBTB7A"]
        tasks = [pool.apply_async(ko_simulation, args=(oracle, tf)) for tf in tf_of_interest]
        results = []
        for task in tasks:
            try:
                result = task.get()  # Specify timeout according to needs
                results.append(result)
                print(f"Task completed successfully for {task}")
            except Exception as e:
                print(f"Task failed: {e}")
                results.append(None)



