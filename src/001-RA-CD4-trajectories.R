#\\\
#\\\
# RUN IN THIS SCRIPT WITH SLURM
#\\\
#\\\

'logs' %>% {if(!dir.exists(.)) dir.create(.)}
if(FALSE&FALSE){

slurm_cmd <- "#!/bin/bash

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
dir=\"/projects/canderson2@xsede.org/trajectory-alignment\"
cd \"$dir\" || exit 1

# --- Load Anaconda and activate environment ---
module purge
module load anaconda
conda activate r4.4_env

# --- Run R script ---
R_SCRIPT=\"$dir/src/001-RA-CD4-trajectories.R\"
cat \"$R_SCRIPT\"
Rscript --vanilla \"$R_SCRIPT\"
"
  # write to file (text first, filename second)
  writeLines(slurm_cmd, "src/001-slurm.sh")
  
  # RUN
  system("sbatch src/001-slurm.sh")
}

#\\\
#\\\
# Set Up
#\\\
#\\\
rm(list = ls()); gc()

pacman::p_load(readr, dplyr, snakecase, magrittr, ggplot2, patchwork, slingshot,
    SingleCellExperiment, TrajectoryUtils, scran, scater, scuttle, SummarizedExperiment, harmony,
    tidyr, tibble, stringr, zeallot,
    install = F)

wd <- readLines('.wd.txt')
setwd(wd)
getwd()

avail <- as.numeric(system('nproc', intern = T) )


#\\\\
#\\\\
# load data
#\\\\
#\\\\

sce <- qs::qread("/projects/canderson2@xsede.org/zhang-lab/cite-seq/analysis-versions/version002/cd4/processed-data/sce002.qs", 
    nthreads = 3)


#\\\\
#\\\\
#\\\\
# Partition data
#\\\\
#\\\\
#+ Subset patients and features to independent datasets

# sample participants maintaining cell-type, visit, and converterstatus
nProp <- 0.5 
set.seed(12039)  
samp_cols <- data.frame(colData(sce)) %>%  
  select(participant_id, VisitType_, new_cluster) %>%  
  rownames_to_column("samp") %>%  
  group_by(participant_id, VisitType_,new_cluster) %>%  
  sample_frac(nProp, replace = FALSE) 

sce_x <- sce[, colnames(sce) %in% samp_cols$samp] 
sce_y <- sce[, !colnames(sce) %in% samp_cols$samp] 

ncol(sce)==(ncol(sce_x)+ncol(sce_y))

# sample genes
pProp = .5
samp_rows  <-  sample(rownames(sce), floor(pProp*nrow(sce)), replace = FALSE)

sce_x <- sce[ rownames(sce) %in% samp_rows , ] 
sce_y <- sce[!rownames(sce) %in% samp_rows , ] 




#\\\\
#\\\\
# Re Run dimension reduction on each 
#\\\\
#\\\\

getPCs <- function(an_sce,Ngenes = NA, Pgenes = NA){
    message("Finding Most Variable Genes")
    # find most variable genes 
    ## model the variance
    dec <- modelGeneVar(an_sce) 
    fit <- metadata(dec)

    ## grab top proportion/N
    if(!is.na(Ngenes)){
        tops <- getTopHVGs(dec, n = Ngenes)
    }else if(!is.na(Pgenes)){
        tops <- getTopHVGs(dec, prop = Pgenes)
    } else{
        warning("Please provide a Ngenes or Pgenes to subset genes with")
    }
    message("Done\n")

    message("Calculating PCs")
    # calculate PCA on most variable genes
    ## Gene PCA for top 20 PCs (uses Irlba which has random start so use seed)
    pca_pars <- list(
        exprs_values = "logcounts",
        scale = TRUE, 
        name = "PCA",
        ncomponents = 50,
        subset_row = tops # subset for these genes to run in pca
    )

    set.seed(1020)
    an_sce <- do.call(scater::runPCA, append(pca_pars, list(x = an_sce)))
    message("Done\n")

    message("Running Harmony Batch Correction")
    # harmony batch correct PCs
    if("HARMONY" %in% names(reducedDims(an_sce))) reducedDims(an_sce)$HARMONY <- NULL

    ## run harmony
    harmony_pars <- list(
        reduction.use = "PCA",
        group.by.vars = c("participant_id"),
        max_iter = 10, # aka harmony.iter
        plot_convergence = TRUE,
        .options = harmony_options(max.iter.cluster = 30,
                                    epsilon.cluster = -Inf,
                                    epsilon.harmony = -Inf),
        ncores = avail-1,
        verbose = TRUE
    )

    set.seed(2303)
    an_sce <- do.call(RunHarmony, append(list(object = an_sce), harmony_pars))
    message("Done\n")

    return(list(an_sce, tops))
}

c(sce_x, top_genes_x) %<-% getPCs(sce_x, Ngenes = 2000)
c(sce_y, top_genes_y) %<-% getPCs(sce_y, Ngenes = 2000)

metadata(sce_x)$metadata$top_genes = top_genes_x
metadata(sce_y)$metadata$top_genes = top_genes_y

#\\\\
#\\\\
# Save
#\\\\
#\\\\

qsave(list(sce_x,sce_y), "processed-data/sce-x-and-y.qs", preset = 'high')