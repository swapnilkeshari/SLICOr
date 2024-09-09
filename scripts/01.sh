#!/bin/bash
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 6-00:00:00 # Runtime in D-HH:MM
#SBATCH --job-name=01    # create a short name for your job
#SBATCH --nodes=1                # node count
##SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=swk25@pitt.edu
#SBATCH --cpus-per-task=16

module purge
module load gcc/12.2.0
module load r/4.3.0
pwd; hostname; date
Rscript 01_atacdata_analysis_with_cicero_and_monocle3.R
crc-job-stats