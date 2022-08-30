#!/bin/bash
#SBATCH --job-name="taxa_age_tess_run100"
#SBATCH --time=200:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # CHANGE THIS to processor core(s) per node
#SBATCH --mail-user="erenada@uri.edu" #CHANGE THIS to your user email address
#SBATCH --mail-type=ALL
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err
#SBATCH --exclusive

module purge

module load R/4.2.0-foss-2021b

Rscript taxa_age_tess_loop_chr100.R
