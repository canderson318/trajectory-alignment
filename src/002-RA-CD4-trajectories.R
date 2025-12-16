
rm(list = ls()); gc()

pacman::p_load(readr, dplyr, snakecase, magrittr, ggplot2, patchwork, slingshot, SingleCellExperiment, TrajectoryUtils, scater, SummarizedExperiment, tidyr, tibble, stringr, zellkonverter)


setwd("~/Documents/personal/R stuff/traj-alignment/")


# # download data
# if(FALSE)"curl https://covid19.cog.sanger.ac.uk/submissions/release2/meyer_nikolic_covid_pbmc_raw.h5ad --output raw-data/meyer_nikolic_covid_pbmc_raw.h5ad" %>% system()


# load data
sce <- readH5AD("raw-data/meyer_nikolic_covid_pbmc_raw.h5ad")

# Check structure
names(assays(sce)) <- "logcounts"
sce <- runPCA(sce, ncomponents = 20)


plotPCA(sce, color_by = )
