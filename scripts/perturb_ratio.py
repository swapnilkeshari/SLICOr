import numpy as np, pandas as pd, scanpy as sc, anndata as ad, matplotlib.pyplot as plt, os, logging, seaborn as sns, pickle
from scipy.stats import median_abs_deviation, hypergeom
import multiprocessing as mp, celloracle as co

# oracle = co.load_hdf5("../out_data/files/oracle_fitted.celloracle.oracle")

def perturb_in_out(goi):
    result = {goi: []}
    print(f"Calculating ratio for {goi}", flush=True)
    oracle = co.load_hdf5(os.path.join(wd, '../out_data/intermediate_data', f"{goi}_perturbation.celloracle.oracle"))
    color_dict = {str(k): v for k, v in zip(set(oracle.adata.obs['leiden'].astype(int)), ['lightskyblue', 'dodgerblue', 'mediumorchid', 'limegreen', 'darkblue', 'darkgray', 'green', 'firebrick', 'sandybrown', 'lightcoral', 'teal', 'gold'])}
    oracle.colorandum = np.array([color_dict[i] for i in oracle.adata.obs['leiden']])    
    cluster_list = list(set(oracle.adata.obs.leiden))
    # Annotation 1
    annotation = {"Naive":[5],
                "GC": [3,6],
                "PB":[7,9],  
                "ABC":[1],
                "Day 1 cells" :[0],   
                "Day3_4":[2,4,8],
                "Undefined":[10,11]}

    # Change dictionary format for annotation assignment
    anno_reverse = {}
    for i in cluster_list:
        for k in annotation:
            if int(i) in annotation[k]:
                anno_reverse[i] = k
    oracle.adata.obs["cell_type"] = [anno_reverse[i] for i in oracle.adata.obs.leiden]


    gc_mean = np.mean(oracle.adata.X[(oracle.adata.obs['cell_type'] == 'GC'),(oracle.adata.var_names == goi)].A)
    pb_mean = np.mean(oracle.adata.X[(oracle.adata.obs['cell_type'] == 'PB'),(oracle.adata.var_names == goi)].A)
    if gc_mean > pb_mean:
        TF_type = 'GC'
    else:
        TF_type = 'PB'
    
    totalGC_counts = oracle.adata.obs[(oracle.adata.obs['cell_type'] == 'GC')].shape[0]
    totalPB_counts = oracle.adata.obs[(oracle.adata.obs['cell_type'] == 'PB')].shape[0]
    total34_counts = oracle.adata.obs[(oracle.adata.obs['cell_type'] == 'Day3_4')].shape[0]

    transition_matrix = oracle.transition_prob
    # transition_matrix = oracle.transition_prob - (oracle.embedding_knn.A/ oracle.embedding_knn.sum(1).A.T)
    transition_matrix = abs(transition_matrix[oracle.adata.obs['cell_type'].isin(['GC', 'PB', 'Day3_4']),:])
    row_labels = oracle.adata.obs[oracle.adata.obs['cell_type'].isin(['GC', 'PB', 'Day3_4'])]['cell_type'].values
    col_labels = oracle.adata.obs['cell_type'].values
    tm_df = pd.DataFrame(transition_matrix, index=row_labels, columns=col_labels)
    # Based on Max probability
    mp = tm_df.idxmax(axis=1).reset_index()
    mp = mp[mp['index']!=mp[0]].groupby('index').agg(lambda x: x.value_counts().to_dict())
    gc_in = 0
    try:
        gc_in += mp.loc['Day3_4'][0]['GC']
    except KeyError:
        pass
    
    try:
        gc_in += mp.loc['PB'][0]['GC']
    except KeyError:
        pass

    pb_in = 0
    try:
        pb_in += mp.loc['Day3_4'][0]['PB']
    except KeyError:
        pass

    try:
        pb_in += mp.loc['GC'][0]['PB']
    except KeyError:
        pass
    
    try:
        gc_out = 0 + sum(list(mp.loc['GC'][0].values()))
    except KeyError:
        gc_out = 0
    
    try:
        pb_out = 0 + sum(list(mp.loc['PB'][0].values()))
    except KeyError:
        pb_out = 0
    # gc_out = 0 + sum(list(mp.loc['GC'][0].values()))
    # pb_out = 0 + sum(list(mp.loc['PB'][0].values()))
    
    finalGC_counts = totalGC_counts - gc_out + gc_in
    finalPB_counts = totalPB_counts - pb_out + pb_in
    result[goi] = [TF_type, totalGC_counts, totalPB_counts, total34_counts, finalGC_counts, finalPB_counts,gc_in, pb_in, gc_out, pb_out]

    U, V =  oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'Naive', 0], oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'Naive', 1]
    magnitudes_Naive = np.sqrt(V**2 + U**2)
    U, V =  oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'Naive', 0], oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'Naive', 1]
    magnitudes_Naive_random = np.sqrt(V**2 + U**2)
    if TF_type == 'PB':
        U, V =  oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'PB', 0], oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'PB', 1]
        magnitudes_PB = np.sqrt(V**2 + U**2)
        U, V =  oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'PB', 0], oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'PB', 1]
        magnitudes_PB_random = np.sqrt(V**2 + U**2)
        data = pd.DataFrame({
            'Magnitude': np.concatenate([magnitudes_Naive, magnitudes_Naive_random, magnitudes_PB, magnitudes_PB_random]),
            'Type': ['Actual']*len(magnitudes_Naive) + ['Random']*len(magnitudes_Naive_random)+ ['Actual']*len(magnitudes_PB) + ['Random']*len(magnitudes_PB_random),
            'cell_type': ['Naive']*len(magnitudes_Naive) + ['Naive']*len(magnitudes_Naive_random)+['PB']*len(magnitudes_PB) + ['PB']*len(magnitudes_PB_random)
        })
    else:
        U, V =  oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'GC', 0], oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'GC', 1]
        magnitudes_GC = np.sqrt(V**2 + U**2)
        U, V =  oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'GC', 0], oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'GC', 1]
        magnitudes_GC_random = np.sqrt(V**2 + U**2)
        data = pd.DataFrame({
            'Magnitude': np.concatenate([magnitudes_Naive, magnitudes_Naive_random, magnitudes_GC, magnitudes_GC_random]),
            'Type': ['Actual']*len(magnitudes_Naive) + ['Random']*len(magnitudes_Naive_random)+ ['Actual']*len(magnitudes_GC) + ['Random']*len(magnitudes_GC_random),
            'cell_type': ['Naive']*len(magnitudes_Naive) + ['Naive']*len(magnitudes_Naive_random)+['GC']*len(magnitudes_GC) + ['GC']*len(magnitudes_GC_random)
        })
    # Plot the actual data on the left half
    sns.violinplot(data=data, x = 'cell_type', y='Magnitude', hue= "Type",split=True, fill=False,inner="quart", linewidth=1)
    plt.savefig(os.path.join(wd, '../out_data/perturbed_data', f'{goi}_{TF_type}_magnitude_violin.png'))
    plt.close()

    n_grid = 40; min_mass = 0.04
    oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
    fig, ax = plt.subplots(figsize=[8, 8])
    scale_simulation = 15
    oracle.plot_cluster_whole(ax=ax, s=10)
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
    plt.title(f"{goi} KO simulation")
    plt.savefig(os.path.join(wd, '../out_data/figures', f"{goi}_KO_simulation.svg"), dpi=600, format="svg")
    plt.close()

    return result

def perturb_tp(goi):
    result = {goi: []}
    print(f"Calculating ratio for {goi}", flush=True)
    oracle = co.load_hdf5(os.path.join(wd, '../out_data/intermediate_data', f"{goi}_perturbation.celloracle.oracle"))
    gc_mean = np.mean(oracle.adata.X[(oracle.adata.obs['cell_type'] == 'GC'),(oracle.adata.var_names == goi)].A)
    pb_mean = np.mean(oracle.adata.X[(oracle.adata.obs['cell_type'] == 'PB'),(oracle.adata.var_names == goi)].A)
    if gc_mean > pb_mean:
        TF_type = 'GC'
    else:
        TF_type = 'PB'
    
    totalGC_counts = oracle.adata.obs[(oracle.adata.obs['cell_type'] == 'GC')].shape[0]
    totalPB_counts = oracle.adata.obs[(oracle.adata.obs['cell_type'] == 'PB')].shape[0]
    transition_matrix = oracle.transition_prob
    # transition_matrix = oracle.transition_prob - (oracle.embedding_knn.A/ oracle.embedding_knn.sum(1).A.T)
    transition_matrix = abs(transition_matrix[((oracle.adata.obs['cell_type'] == 'GC') | (oracle.adata.obs['cell_type'] == 'PB')),:])
    row_labels = oracle.adata.obs[((oracle.adata.obs['cell_type'] == 'GC') | (oracle.adata.obs['cell_type'] == 'PB'))]['cell_type'].values
    col_labels = oracle.adata.obs['cell_type'].values
    tm_df = pd.DataFrame(transition_matrix, index=row_labels, columns=col_labels)
    result[goi] = [TF_type, totalGC_counts, totalPB_counts]

    # Based on Max probability
    mp = tm_df.idxmax(axis=1).reset_index()
    # Based on sum of probabilities
    sp = tm_df.T.groupby(tm_df.T.index).sum().T.idxmax(axis=1).reset_index()
    # Based on max counts above threshold
    tp = (tm_df.T>0.1).groupby(tm_df.T.index).sum().T
    tp_rowsum = tp.sum(axis=1)
    final_cell=[]
    for i,start_cell in enumerate(tp.index):
        if tp_rowsum.iloc[i] == 0:
            final_cell.append(start_cell)
        else:
            final_cell.append(tp.iloc[i].idxmax())
    tp_df = pd.DataFrame()
    tp_df['index'] = tp.index
    tp_df[0] = final_cell
    for df in [mp,sp,tp_df]:
        exitGC_counts = sum(df[df['index'] == 'GC']['index'] != df[df['index'] == 'GC'][0])
        exitPB_counts = sum(df[df['index'] == 'PB']['index'] != df[df['index'] == 'PB'][0])
        result[goi].extend([exitGC_counts, exitPB_counts])
    

    U, V =  oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'Naive', 0], oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'Naive', 1]
    magnitudes_Naive = np.sqrt(V**2 + U**2)
    U, V =  oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'Naive', 0], oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'Naive', 1]
    magnitudes_Naive_random = np.sqrt(V**2 + U**2)
    if TF_type == 'PB':
        U, V =  oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'PB', 0], oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'PB', 1]
        magnitudes_PB = np.sqrt(V**2 + U**2)
        U, V =  oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'PB', 0], oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'PB', 1]
        magnitudes_PB_random = np.sqrt(V**2 + U**2)
        data = pd.DataFrame({
            'Magnitude': np.concatenate([magnitudes_Naive, magnitudes_Naive_random, magnitudes_PB, magnitudes_PB_random]),
            'Type': ['Actual']*len(magnitudes_Naive) + ['Random']*len(magnitudes_Naive_random)+ ['Actual']*len(magnitudes_PB) + ['Random']*len(magnitudes_PB_random),
            'cell_type': ['Naive']*len(magnitudes_Naive) + ['Naive']*len(magnitudes_Naive_random)+['PB']*len(magnitudes_PB) + ['PB']*len(magnitudes_PB_random)
        })
    else:
        U, V =  oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'GC', 0], oracle.delta_embedding[oracle.adata.obs['cell_type'] == 'GC', 1]
        magnitudes_GC = np.sqrt(V**2 + U**2)
        U, V =  oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'GC', 0], oracle.delta_embedding_random[oracle.adata.obs['cell_type'] == 'GC', 1]
        magnitudes_GC_random = np.sqrt(V**2 + U**2)
        data = pd.DataFrame({
            'Magnitude': np.concatenate([magnitudes_Naive, magnitudes_Naive_random, magnitudes_GC, magnitudes_GC_random]),
            'Type': ['Actual']*len(magnitudes_Naive) + ['Random']*len(magnitudes_Naive_random)+ ['Actual']*len(magnitudes_GC) + ['Random']*len(magnitudes_GC_random),
            'cell_type': ['Naive']*len(magnitudes_Naive) + ['Naive']*len(magnitudes_Naive_random)+['GC']*len(magnitudes_GC) + ['GC']*len(magnitudes_GC_random)
        })
    # Plot the actual data on the left half
    sns.violinplot(data=data, x = 'cell_type', y='Magnitude', hue= "Type",split=True, fill=False,inner="quart", linewidth=1)
    plt.savefig(os.path.join(wd, '../out_data/perturbed_data', f'{goi}_{TF_type}_magnitude_violin.png'))
    plt.close()
    return result


def classify_TF_perturb(goi):
    print(f"Calculating ratio for {goi}", flush=True)
    oracle = co.load_hdf5(os.path.join(wd, '../out_data/intermediate_data', f"{goi}_perturbation.celloracle.oracle"))
    gc_mean = np.mean(oracle.adata.X[(oracle.adata.obs['cell_type'] == 'GC'),(oracle.adata.var_names == goi)].A)
    pb_mean = np.mean(oracle.adata.X[(oracle.adata.obs['cell_type'] == 'PB'),(oracle.adata.var_names == goi)].A)
    
    if gc_mean > pb_mean:
        TF_type = 'GC'
    else:
        TF_type = 'PB'
    
    delta_embed = pd.DataFrame(oracle.delta_embedding)
    delta_embed.columns = ['delta_x', 'delta_y']
    delta_embed['magnitude'] = delta_embed.apply(lambda x: np.sqrt(x['delta_x']**2 + x['delta_y']**2), axis=1)
    delta_embed['cell_type'] = oracle.adata.obs['cell_type'].values
    delta_embed['angles'] = delta_embed.apply(lambda x: np.degrees(np.arctan2(x['delta_y'], x['delta_x'])), axis=1)
    # angles_GC, mag_GC = delta_embed[delta_embed['cell_type'] == 'GC']['angles'].values, delta_embed[delta_embed['cell_type'] == 'GC']['magnitude'].values
    # angles_PB, mag_PB = delta_embed[delta_embed['cell_type'] == 'PB']['angles'].values, delta_embed[delta_embed['cell_type'] == 'PB']['magnitude'].values
    median_naive = np.median(delta_embed[delta_embed['cell_type'] == 'Naive']['magnitude'].values)
    totalGC_counts = delta_embed[(delta_embed['cell_type'] == 'GC')].shape[0]
    totalPB_counts = delta_embed[(delta_embed['cell_type'] == 'PB')].shape[0]
    if TF_type == 'GC':
        exitGC_counts = delta_embed[(delta_embed['cell_type'] == 'GC') & ((delta_embed['angles']<0) & (delta_embed['angles']> -90)) & (delta_embed['magnitude'] > max(median_naive, 0.3))].shape[0]
        exitPB_counts = delta_embed[(delta_embed['cell_type'] == 'PB') & ~((delta_embed['angles']<0) & (delta_embed['angles']> -90)) & (delta_embed['magnitude'] > max(median_naive, 0.3))].shape[0]
    else:
        exitGC_counts = delta_embed[(delta_embed['cell_type'] == 'GC') & ~(delta_embed['angles'] > 90) & (delta_embed['magnitude'] > max(median_naive, 0.3))].shape[0]
        exitPB_counts = delta_embed[(delta_embed['cell_type'] == 'PB') & (delta_embed['angles'] > 90) & (delta_embed['magnitude'] > max(median_naive, 0.3))].shape[0]

    return goi, TF_type, exitGC_counts,totalGC_counts, exitPB_counts, totalPB_counts,((totalGC_counts-exitGC_counts )/ (totalPB_counts - exitPB_counts)),(totalGC_counts/totalPB_counts)

def ratio_enrich(goi):
    try:
        print(f"Calculating ratio for {goi}", flush=True)
        oracle = co.load_hdf5(os.path.join(wd, '../out_data/intermediate_data', f"{goi}_perturbation.celloracle.oracle"))
        delta_embed = pd.DataFrame(oracle.delta_embedding)
        delta_embed.columns = ['delta_x', 'delta_y']
        delta_embed['magnitude'] = pd.DataFrame(np.linalg.norm(oracle.delta_embedding, ord=2, axis=1))
        delta_embed['cell_type'] = oracle.adata.obs['cell_type'].values
        delta_embed['slope'] = (delta_embed['delta_y']/delta_embed['delta_x'])#.apply(lambda x: math.degrees(math.atan(x)))
        proGC_counts = delta_embed[(delta_embed['cell_type'] == 'PB') & (delta_embed['slope'] < 0) & (delta_embed['delta_x'] < 0)].shape[0]
        proPB_counts = delta_embed[(delta_embed['cell_type'] == 'GC') & (delta_embed['slope'] <0) & (delta_embed['delta_y'] < 0)].shape[0]
        median_naive = np.median(delta_embed['magnitude'][delta_embed['cell_type'] == 'Naive'])
        exitGC_counts = delta_embed[(delta_embed['cell_type'] == 'GC') & (delta_embed['magnitude'] > median_naive)].shape[0]
        exitPB_counts = delta_embed[(delta_embed['cell_type'] == 'PB') & (delta_embed['magnitude'] > median_naive)].shape[0]
        exitMBC_counts = delta_embed[(delta_embed['cell_type'] == 'MBC from GC') & (delta_embed['magnitude'] > median_naive)].shape[0]
        totalGC_counts = len(delta_embed[(delta_embed['cell_type'] == 'GC')])
        totalPB_counts = len(delta_embed[(delta_embed['cell_type'] == 'PB')])
        totalMBC_counts = len(delta_embed[(delta_embed['cell_type'] == 'MBC from GC')])
        return (goi, proGC_counts, proPB_counts, exitGC_counts, exitPB_counts, exitMBC_counts, totalGC_counts, totalPB_counts, totalMBC_counts)
    except Exception as e:
        print(f'Process {mp.current_process().name} failed for gene {goi} with error {e}', flush=True)
        return (goi, -1, -1, -1, -1, -1, -1, -1, -1)

if __name__ == '__main__':
    plt.rcParams["savefig.dpi"] = 300
    wd = '/ix/djishnu/Swapnil/CellOracle/BCell/primaryBCell/notebooks'
    out_path = os.path.join(wd, '../out_data')
    os.makedirs(f"{out_path}/figures", exist_ok=True)
    os.makedirs(f"{out_path}/intermediate_data", exist_ok=True)
    os.makedirs(f"{out_path}/perturbed_data", exist_ok=True)
    sc.settings.figdir = f"{out_path}/figures"
    # ko_simulation(oracle, "BCL6")
    with mp.Pool(processes=16) as pool:
        tf_of_interest = ["BATF", "IRF4", "IRF8", "PRDM1", "SPIB"]
        # tf_of_interest = ["BATF3","BCL6","EGR1","FOS","IKZF1","IRF1","JUNB","JUND","MEF2A","MEF2C","MYB","NFATC1","NFATC2","NFIL3","NFKB1","NFKB2","PAX5","RUNX1","SP3","SPI1","STAT1","STAT5A","TCF12","VDR","XBP1","ZBTB7A", "BATF", "IRF4", "IRF8", "PRDM1", "SPIB"]

        results = pool.map(perturb_in_out, tf_of_interest)
        # pd.DataFrame(results, columns = ['goi', 'TF_type', 'exitGC_counts', 'totalGC_counts', 'exitPB_counts', 'totalPB_counts', 'Final_Ratio', 'Initial_Ratio']).to_csv(os.path.join(out_path, 'intermediate_data', 'ratio_enrichment.csv'), index=False)
        with open(os.path.join(out_path, 'intermediate_data', 'ratio_enrichment.pkl'), 'wb') as f:
            pickle.dump(results, f)
        # for task in tasks:
        #     try:
        #         result = task.get()  # Specify timeout according to needs
        #         print(f"Task completed: {task}")
        #     except Exception as e:
        #         print(f"Task failed: {e}")