#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00:00 # Runtime in D-HH:MM
#SBATCH --job-name=05_multiome    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=4    # total number of tasks across all nodes
##SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=swk25@pitt.edu
#SBATCH --cpus-per-task=4        # cpu-cores per task
#SBATCH --mem=800G
#SBATCH --cluster=smp
#SBATCH --partition=high-mem
#SBATCH --constraint=intel,broadwell

module purge
module load gcc/8.2.0
module load python/anaconda3.8-2020.11
# module load bedtools/2.30.0
source activate celloracle_env

python 05_KO_1.py
echo "bashdone"
crc-job-stats