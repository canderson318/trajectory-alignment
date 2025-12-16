#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --ntasks=10
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --qos=normal

#SBATCH --job-name=001-job

#SBATCH --output=/projects/canderson2@xsede.org/trajectory-alignment/logs/001-slurm-%j.out
#SBATCH --error =/projects/canderson2@xsede.org/trajectory-alignment/logs/001-slurm-%j.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=christian.anderson@cuanschutz.edu


# --- Set working directory ---
dir="/projects/canderson2@xsede.org/trajectory-alignment"
cd "$dir" || exit 1

# --- Load Anaconda and activate environment ---
module purge
module load anaconda
conda activate r4.4_env

# --- Run R script ---
R_SCRIPT="$dir/src/001-RA-CD4-trajectories.R"
cat "$R_SCRIPT"
Rscript --vanilla "$R_SCRIPT"

