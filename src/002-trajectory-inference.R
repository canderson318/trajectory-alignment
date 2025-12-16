#\\\
#\\\
# RUN IN THIS SCRIPT WITH SLURM
#\\\
#\\\

# system("sbatch src/002-slurm.sh")


#\\\
#\\\
# Set Up
#\\\
#\\\
rm(list = ls()); gc()

pacman::p_load(readr, dplyr, snakecase, magrittr, ggplot2, patchwork, slingshot,
    SingleCellExperiment, TrajectoryUtils, scran, scater, scuttle, SummarizedExperiment, harmony,
    tidyr, tibble, stringr, zeallot, parallel,doParallel, myPackage
    install = F)

wd <- readLines('~/.wd.txt')
setwd(wd)
getwd()

avail <- as.numeric(system('nproc', intern = T) )

# Detect and use # of cores 
cores <- avail - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

c(sce_x, sce_y)%<-% qs::qread("processed-data/sce-x-and-y.qs", nthreads = 3)



# \\\
# \\\
# Run Slingshot on each dataset
# \\\
# \\\
calc.slingshot <- function(X, labels){
  slingshot::slingshot(data = X, 
                        clusterLabels = labels, 
                        omega = FALSE, 
                        use.median = FALSE, 
                        approx_points = 10,
                        start.clus = "Naive CD4"
                        )
}

args = list('x' = list(X = reducedDims(sce_x)$PCA, labels = sce_x$new_cluster)), 
            'y' = list(X = reducedDims(sce_y)$PCA, labels = sce_y$new_cluster))

c(sling_x, sling_y) %<-% mclapply(args, function(z){
    X = z$X
    labels = z$labels
    calc.slingshot(X, labels)
})

qs::qsave(list(sling_x, sling_y), 'processed-data/x-and-y-slingshots.qs')

# # extract slingshot curve coordinates
# get.curve.df <- function(curves){
#   curves %>%
#   SlingshotDataSet() %>%
#   slingCurves()  %>%
#     lapply(function(crv) {
#       d <- as.data.frame(crv$s)

#       d
#     }) %>%
#     do.call(rbind, .) %>%
#     rownames_to_column("lineage") %>%
#     mutate(lineage = str_remove(lineage, "\\.\\d+$"))
# }

# # get curves from slingshot objects
# getCurves(...)

# get.curve.df()