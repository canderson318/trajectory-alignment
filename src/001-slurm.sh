#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --qos=normal

#SBATCH --job-name=001-job

#SBATCH --output=/projects/canderson2@xsede.org/trajectory-alignment/logs/001-slurm-%j.out
#SBATCH --error=/projects/canderson2@xsede.org/trajectory-alignment/logs/001-slurm-%j.err

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
rel_script_pth="src/001-RA-CD4-trajectories.R"
R_SCRIPT="$dir/$rel_script_pth"
cat "$R_SCRIPT"
Rscript --vanilla "$R_SCRIPT"

