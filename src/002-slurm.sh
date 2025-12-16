#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --qos=normal

#SBATCH --job-name=002-job

#SBATCH --output=/projects/canderson2@xsede.org/trajectory-alignment/logs/002-slurm-%j.out
#SBATCH --error=/projects/canderson2@xsede.org/trajectory-alignment/logs/002-slurm-%j.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=christian.anderson@cuanschutz.edu


# --- Set working directory ---
dir="/projects/canderson2@xsede.org/trajectory-alignment"
cd "$dir" || exit 1

# --- Make loggind directory ---
mkdir -p "$dir/logs"


# --- Load Anaconda and activate environment ---
module purge
module load anaconda
conda activate r4.4_env

# --- Run R script ---
rel_script_pth="src/002-trajectory-inference.R"
R_SCRIPT="$dir/$rel_script_pth"
cat "$R_SCRIPT"
Rscript --vanilla "$R_SCRIPT"

