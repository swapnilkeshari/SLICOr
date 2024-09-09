#!/bin/bash
#SBATCH --job-name="03"
#SBATCH --output="/ocean/projects/cis240075p/skeshari/igvf/bcell2/primaryBCell/scripts/03.log"
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH --ntasks-per-node=32
#SBATCH -t 0-01:30:00

module load anaconda3/2022.10
source activate celloracle_env

python /ocean/projects/cis240075p/skeshari/igvf/bcell2/primaryBCell/scripts/03_scanpy_preprocessing.py