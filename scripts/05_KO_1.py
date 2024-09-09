import numpy as np, pandas as pd, scanpy as sc, anndata as ad, matplotlib.pyplot as plt, os, logging, seaborn as sns, pickle
from scipy.stats import median_abs_deviation, hypergeom
import multiprocessing as mp, celloracle as co

# oracle = co.load_hdf5("../out_data/files/oracle_fitted.celloracle.oracle")
def ko_simulation(oracle,goi, cluster_column_name="cell_type"):
    try:
        # Enter perturbation conditions to simulate signal propagation after the perturbation.
        print(f"Simulating perturbation for {goi}")
        oracle.simulate_shift(perturb_condition={goi: 0.0},n_propagation=2) # If kernal crashes increase the memory limit of jupyter notebook
        print(f"Estimating transition probability for {goi}")
        oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1) # Get transition probability
        print(f"Calculating embedding for {goi}")
        oracle.calculate_embedding_shift(sigma_corr=0.05) # Calculate embedding
        print(f'Saving perturbed data for {goi}')
        oracle.to_hdf5(os.path.join(wd, '../out_data/intermediate_data',f'{goi}_perturbation.celloracle.oracle'))
    except Exception as e:
        print(f'Process {mp.current_process().name} failed for gene {goi} with error {e}', flush=True)
    return None

if __name__ == '__main__':
    plt.rcParams["savefig.dpi"] = 300
    wd = '/ix/djishnu/Swapnil/CellOracle/BCell/primaryBCell/notebooks'
    out_path = os.path.join(wd, '../out_data')
    os.makedirs(f"{out_path}/figures", exist_ok=True)
    os.makedirs(f"{out_path}/intermediate_data", exist_ok=True)
    os.makedirs(f"{out_path}/pertubed_data", exist_ok=True)
    sc.settings.figdir = f"{out_path}/figures"

    oracle = co.load_hdf5(os.path.join(wd, '../out_data/intermediate_data', "oracle_fitted.celloracle.oracle"))

    # ko_simulation(oracle, "BCL6")
    with mp.Pool(processes=4) as pool:
        # tf_of_interest = ["BATF", "IRF4", "IRF8", "PRDM1", "SPIB"]
        tf_of_interest = ["BATF3","BCL6","EGR1","FOS","IKZF1","IRF1","JUNB","JUND","MEF2A","MEF2C","MYB","NFATC1","NTC1","NFATC2","NFIL3","NFKB1","NFKB2","PAX5","","RUNX1","SP3","SPI1","STAT1","STAT5A","TCF12","VDR","XBP1","ZBTB7A"]

        tasks = [pool.apply_async(ko_simulation, args=(oracle, tf)) for tf in tf_of_interest]
        for task in tasks:
            try:
                result = task.get()  # Specify timeout according to needs
                print(f"Task completed: {task}")
            except Exception as e:
                print(f"Task failed: {e}")