#\\\
#\\\
# Set up
#\\\
#\\\
rm(list = ls()); gc()

pacman::p_load(readr, dplyr, snakecase, magrittr, ggplot2, patchwork, slingshot, SingleCellExperiment, TrajectoryUtils, scater, SummarizedExperiment, tidyr, tibble, stringr, install = F)

library(myPackage)

# wd <- "/projects/canderson2@xsede.org/trajectory-alignment"
# wd <- "/Users/canderson/trajectory-alignment"
# writeLines(wd, '~/.wd.txt')

wd <- readLines('~/.wd.txt')
setwd(wd)
wd

tmp_plot  <-  function( path = sprintf("%s/results/temp-plot.png", wd), ...){
  temp_plot(path = path, ...)
}

#\\\
#\\\
# Download Data
#\\\
#\\\

# make raw-data directory 
"raw-data" %>% 
  {if(!dir.exists(.)) dir.create(.)}

if(!file.exists('raw-data/wine.data')){
  system("curl https://archive.ics.uci.edu/static/public/109/wine.zip --output raw-data/wine.zip")
  system("ls")
  system("unzip raw-data/wine.zip -d raw-data")
  system('rm raw-data/wine.zip')
}

system("cat raw-data/wine.data")
system("cat raw-data/wine.names")

#\\\
#\\\
# Load Data
#\\\
#\\\

nms <- snakecase::to_snake_case(c("class","Alcohol","Malic acid","Ash","Alcalinity of ash","Magnesium","Total phenols","Flavanoids","Nonflavanoid phenols","Proanthocyanins","Color intensity","Hue","OD280/OD315 of diluted wines","Proline"))
wine <- read_csv("raw-data/wine.data", col_names = nms)
wine$class <- c("A", "B", "C")[wine$class]


#\\\
#\\\
# Scale Data
#\\\
#\\\

scaled <- scale(wine[,-1])
scaled <- cbind(wine[,1], scaled)

# make data object to store all objects
obj <- list(data = list(orig = wine, scaled = scaled), reducedDims = list(pcs = NA))

#\\\
#\\\
# Calculate PCs on whole dataset
#\\\
#\\\
calc.pcs <- function(X){
  z <- prcomp(X)
  txtplot::txtplot(z$sdev^2 %>% {(./sum(.))*100})
  return(z$x)
}

# plot both pcs
x <-calc.pcs(obj$data$scaled[, -1])

ab_pcs <- ggplot(x, aes(PC1, PC2, color = wine[,1, drop = TRUE] %>% factor()))+
  stat_ellipse(type = "norm", alpha = .5, level = .95)+
  geom_point(size = 3, alpha = .5)+
  theme_bw()+
  labs(color = "class")+
  scale_color_hue()+
  theme(aspect.ratio = 1)+
  ggtitle("Dataset AB")


# \\\
# \\\
# Partition Data
# \\\
# \\\

set.seed(18)
# # subset columnwise
# a_cols <- sample(2:14, floor(14/2))
# b_cols <- setdiff(2:14, a_cols)

# subset columns with some shared features
a_cols <- 2:10
b_cols <- 6:14

b_rows <- a_rows <- 1:nrow(wine)

# subsample rowwise
a_rows <- sample(1:nrow(wine), floor(nrow(wine)/2))
b_rows <- setdiff(1:nrow(wine), a_rows)

# bootstrap
a_rows <- sample(a_rows, 200, T)
b_rows <- sample(b_rows, 200, T)

obj$data$a <- obj$data$scaled[a_rows,a_cols]
obj$data$b <- obj$data$scaled[b_rows,b_cols]

obj$data$idx <- list(a = list(rows = a_rows, columns = a_cols), 
                     b = list(rows = b_rows, columns = b_cols)
                     )

# \\\
# \\\
# Calculate Reduced Dims on Partitioned Datasets
# \\\
# \\\
obj$reducedDims$pcs <- list(a = calc.pcs(obj$data$a),
                            b = calc.pcs(obj$data$b )
                            )

# \\\
# \\\
# Run Slingshot on each dataset
# \\\
# \\\
calc.slingshot <- function(X, groups, components= NA){
  if(is.na(components)) components = ncol(X)
  slingshot::slingshot(data = X[,1:components], 
                        clusterLabels = groups, 
                        omega = FALSE, 
                        use.median = FALSE, 
                        approx_points = 10,
                        start.clus = "A"
                        )
}

components = NA # number of pcs to use, otherwise uses all
obj$traj <- list(a = calc.slingshot(obj$reducedDims$pcs$a, obj$data$orig$class[a_rows], components =components),
                 b =  calc.slingshot(obj$reducedDims$pcs$b, obj$data$orig$class[b_rows], components =components)
                 )

# get curves from slingshot objects
obj$traj$curves$a <- getCurves(obj$traj$a)
obj$traj$curves$b <- getCurves(obj$traj$b)


# extract slingshot curve coordinates
get.curve.df <- function(curves){
  curves %>%
  SlingshotDataSet() %>%
  slingCurves()  %>%
    lapply(function(crv) {
      d <- as.data.frame(crv$s)

      d
    }) %>%
    do.call(rbind, .) %>%
    rownames_to_column("lineage") %>%
    mutate(lineage = str_remove(lineage, "\\.\\d+$"))
}

obj$traj$curve_dfs <- list(a = get.curve.df(obj$traj$curves$a),
                           b = get.curve.df(obj$traj$curves$b)
                           )

# \\\
# \\\
# Plot Trajectories
# \\\
# \\\
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


xa <- obj$reducedDims$pcs$a%>%
  bind_cols(class= obj$data$scaled[a_rows,1]) %>% 
  data.frame() %>%
  distinct()

xb <- obj$reducedDims$pcs$b%>%
  bind_cols(class= obj$data$scaled[b_rows,1]) %>% 
  data.frame() %>%
  distinct()

a_traj <- ggplot(xa, aes(PC1, PC2))+
  geom_point(aes(color =class %>% factor()), 
    size = 1, alpha = .6)+
  theme_bw()+
  geom_path(data = obj$traj$curve_dfs$a, 
    aes(PC1, PC2, color = lineage),
    color = gg_color_hue(4)[3], linewidth = 1, alpha = .9)+
  labs(color = "class")+
  scale_color_hue()+
  theme(aspect.ratio = 1)+
  ggtitle("Trajectories X and Y")

b_traj <- ggplot(xb, aes(PC1, PC2))+
  geom_point(aes(color =class %>% factor()), size = 1, alpha = .6)+
  geom_path(data = obj$traj$curve_dfs$b, 
    aes(PC1, PC2, color = lineage), 
    color = gg_color_hue(4)[1], linewidth = 1, alpha = .9)+
  labs(color = "class")+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_color_hue()

  
(a_traj/b_traj) %>% tmp_plot( plot = .)

# \\\
# \\\
# Align trajectories with Procrustes
# \\\
# \\\

cols <- seq(2,  do.call(min, lapply(obj$traj$curve_dfs, ncol)) ) # exclude 'lineage' at [1]
# cols <- cols[1:2]

X <- obj$traj$curve_dfs$a[,cols] %>% as.matrix() 
Y <- obj$traj$curve_dfs$b[,cols] %>% as.matrix()

# # Vegan implementation
# procrust <- vegan::procrustes(X = X, Y = Y)
# Y_aligned <- procrust$Yrot

# obj$procrust <- procrust
# (proc_test <- vegan::protest(X, Y))

my_procrustes <- function(X,Y){
  
}

#\\\\
#\\\\
# VISUALIZE Alignment
#\\\\
#\\\\

X_d <- as.data.frame(X)
Y_matched <- as.data.frame(Y_aligned )

colnames(Y_matched) <- paste0("aligned_", cols-1)


b_to_a_traj_points <- ggplot()+

  geom_point(data = obj$reducedDims$pcs$a %>% bind_cols(class_a = obj$data$orig$class[a_rows]) ,aes(PC1, PC2, color = class_a), alpha = .3)+
  
  ggnewscale::new_scale_color()+

  geom_point(data = obj$reducedDims$pcs$b %>% bind_cols(class_b = obj$data$orig$class[b_rows]) ,aes(PC1, PC2, color = class_b), alpha = .3)+
  
  ggnewscale::new_scale_color()+

  geom_path(data =X_d %>% mutate(lineage = "reference"),      aes(PC1, PC2, color = lineage),
            linewidth = 2, alpha = .5)+
  
  geom_path(data = Y_matched %>% mutate(lineage = "aligned"), aes(aligned_1, aligned_2, color = lineage),
            linewidth = 2, alpha = .5)+

  scale_color_hue()+
  ggtitle("Trajectory Y aligned to X with points")+
  theme_bw()

b_to_a_traj_points %>% tmp_plot(plot = .)



b_to_a_traj <- ggplot()+
  
  geom_point(data = obj$reducedDims$pcs$a %>%
               data.frame() %>% mutate(lineage  = "reference") %>% distinct(),
             aes(PC1, PC2, color = lineage),
             alpha = .7)+
  
  geom_point(data = obj$reducedDims$pcs$b %*% procrust$rotation %>%
               data.frame() %>% mutate(lineage  = "aligned") %>% distinct(),
             aes(X1,X2, color = lineage),
             alpha = .7)+
  
  geom_path(data = X_d %>% mutate(lineage = "reference"),      
            aes(PC1, PC2, color = lineage),
            linewidth = 2, alpha = .9)+
  
  geom_path(data = Y_matched %>% mutate(lineage = "aligned"), 
            aes(aligned_1, aligned_2, color = lineage),
            linewidth = 2, alpha = .9)+
  
  scale_color_hue()+
  ggtitle("Trajectory Y aligned to X")+
  labs(color = "")+
  guides(color = guide_legend(
    override.aes = list(
      shape = c(16, 16),   # points
      linetype = c("solid", "solid"),  # lines
      size = c(3, 3)       # larger point legend keys
    )
  )) +
  theme_bw()
b_to_a_traj%>% tmp_plot(plot = .)

p <- (a_traj/b_traj )|(b_to_a_traj+theme(aspect.ratio = 1))

p%>% tmp_plot(plot = .)


