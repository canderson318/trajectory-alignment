#\\
#\\
# SETUP
#\\
#\\

rm(list = ls()); gc()

setwd("~/Documents/personal/R stuff/traj-alignment/")

pacman::p_load(readr, dplyr, snakecase, magrittr, ggplot2, patchwork, slingshot, SingleCellExperiment, TrajectoryUtils, scater, SummarizedExperiment, tidyr, tibble, stringr, myPackage, zeallot)

options(repr.plot.width = 6, repr.plot.height = 6)  # for VSCode/Jupyter envs
par(mfrow = c(1,1))
# or
# dev.new(width = 6, height = 6)
# or
# png("plot.png", width = 6, height = 6, units = "in", res = 300)


#\\
#\\
# LOAD DATA
#\\
#\\

obj <- qs::qread("processed-data/obj001.qs")

# rapply(obj, names, how = "list")

# # two semi dependent matrices with unmatched samples 
# list(X = obj$data$idx$a$rows, Y = obj$data$idx$b$rows) %>% setdiff_intersect()
# list(X = obj$data$idx$a$columns, Y = obj$data$idx$b$columns) %>% setdiff_intersect()

# # making them have unique obs  (bootstrapped in last script)
# X <- obj$data$orig[unique(obj$data$idx$a$rows),obj$data$idx$a$columns ]
# Y <- obj$data$orig[unique(obj$data$idx$b$rows),obj$data$idx$b$columns ]

# two unmatched datasets with same features
set.seed(1030)
X_rows <- sample( 1:nrow(obj$data$orig), floor(nrow(obj$data$orig)/2))
Y_rows <- setdiff(1:nrow(obj$data$orig), X_rows)

X <- obj$data$orig[X_rows,-1] %>% 
  as.matrix()

X_groups <- obj$data$orig$class[X_rows]
Y_groups <- obj$data$orig$class[Y_rows]

Y <- obj$data$orig[Y_rows, -1] %>% 
  as.matrix()

rownames(X) <- X_rows
rownames(Y) <- Y_rows

# # reorder cols
# X <- X[, sample(seq(ncol(X)))]
# Y <- Y[, sample(seq(ncol(Y)))]

# center and scale using X mean and sd
X_mean <- colMeans(X)
X_sd <- colSds(as.matrix(X))

Y_mean <- colMeans(Y)
Y_sd <- colSds(as.matrix(Y))

scldX <- X %>% 
  sweep(., 2, X_mean, "-") %>%
  sweep(., 2, X_sd, "/")
boxplot(scldX, las = 2)

scldY_to_X<- Y %>% 
  sweep(., 2, X_mean, "-") %>% sweep(., 2, pmax(X_sd, .Machine$double.eps), "/")
  
scldY<- Y %>% 
  sweep(., 2, Y_mean, "-") %>% sweep(., 2, Y_sd, "/")
boxplot(scldY, las = 2)



#\\
#\\
# REDUCED DIMENSIONS
#\\
#\\

# SVD on both
## singular values; observation coordinates; axes directions
c(dX_vec, uX, vX) %<-% svd(scldX)
# select top N features
n_comp <- 5
dX_vec <- dX_vec[1:n_comp]
uX <- uX[ , 1:n_comp]
vX <- vX[ , 1:n_comp]
# construct ssingular value matrices
dX <- matrix(0, nrow = length(dX_vec), ncol = length(dX_vec))
diag(dX) <- dX_vec

c(dY_vec, uY, vY) %<-% svd(scldY)
# select top N features
n_comp <- 5
dY_vec <- dY_vec[1:n_comp]
uY <- uY[ , 1:n_comp]
vY <- vY[ , 1:n_comp]
# construct ssingular value matrices
dY <- matrix(0, nrow = length(dY_vec), ncol = length(dY_vec))
diag(dY) <- dY_vec


# sqrt(sum(prcomp(X)$x - uX%*%dX)^2)

# # plot spaces
# par(mfrow = c(3,1))
# plot(uX%*%dX, type = "n", main = "Sample space")
# text(uX%*%dX, labels = rownames(X))

# plot(vX, type = "n", 
#      ylim = c(-1,1), main = "Feature space")
# abline(h = 0, lty =2, col= 'grey50')
# abline(v = 0, lty =2, col= 'grey50')
# text(vX, labels = colnames(X))

# plot(dX_vec  %>% {.^2} %>% {./sum(.)} %>% {.*100}, type = "b", main = "Percent variance explained")
# par(mfrow = c(1,1))

# # project samples into basis space same as sample space in svd
# plot(scldX%*%vX, main = "X•V")
# points(uX%*%dX,main = "U•D", col = "blue", cex = 1.5)
# points(prcomp(scldX)$x,main = "PCA", col = "red", cex = 2)



# \\
# \\
# ALIGN Y basis TO X basis 
# \\
# \\

# scldY_to_X %>% colnames()
# scldX %>% colnames()


X_reduced <- (scldX%*%vX) 
colnames(X_reduced) <- paste0("dim",1:ncol(X_reduced))
dim(X_reduced)

Y_reduced_on_X_v <- (scldY_to_X%*%vX)
colnames(Y_reduced_on_X_v) <- paste0("dim",1:ncol(Y_reduced_on_X_v))
dim(Y_reduced_on_X_v)

# plot(X_reduced, col = "blue")
# points(Y_reduced_on_X_v, col = "red")



### CALCULATE THEIR TRAJECTORIES 

calc.slingshot <- function(X, groups, components= NA){
  if(is.na(components)) components = ncol(X)
  slingshot::slingshot(data = X[,1:components], 
                       clusterLabels = groups, 
                       omega = FALSE, 
                       use.median = FALSE, 
                       approx_points = 50,
                       start.clus = "A"
                       )
}

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


slingX <- calc.slingshot(X_reduced, X_groups)
curveX <- get.curve.df(slingX)

slingY <- calc.slingshot(Y_reduced_on_X_v, Y_groups)
curveY <- get.curve.df(slingY)


Y_scld_to_X_times_X_basis_p <- ggplot()+
  geom_point(data = X_reduced, aes(dim1, dim2), size = 3, color = "salmon")+
  geom_path(data = curveX, aes(dim1,dim2), color = "salmon", linewidth = 3)+
  geom_point(data = Y_reduced_on_X_v, aes(dim1, dim2), size = 3, color = "skyblue")+
  geom_path(data = curveY, aes(dim1,dim2), color = "skyblue", linewidth = 3)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  ggtitle("Y scaled to X * X basis -> curve")

sqrt(mean((as.matrix(curveX[,-1])-as.matrix(curveY[,-1]))^2))
#0.8418248

#\\
#\\
# Align Y BASIS TO X WITH PROCRUSTES THEN CALCULATE TRAJECTORIES
#\\
#\\

vY_aligned_procrust <- procrustes(vX, vY, scale = FALSE)
vY_aligned <- vY_aligned_procrust$Yrot
reducedY_by_rotated_v <- scldY%*%vY_aligned 
colnames(reducedY_by_rotated_v) <- paste0("dim",1:ncol(reducedY_by_rotated_v))

slingY <- calc.slingshot(reducedY_by_rotated_v, Y_groups) 
curveY <- get.curve.df(slingY)

Y_basis_rotated_to_X_basis_p <- ggplot()+
  geom_point(data = X_reduced, aes(dim1, dim2), size = 3, color = "salmon")+
  geom_path(data = curveX, aes(dim1,dim2), color = "salmon", linewidth = 3)+
  geom_point(data = reducedY_by_rotated_v, aes(dim1, dim2), size = 3, color = "skyblue")+
  geom_path(data = curveY , aes(dim1,dim2), color = "skyblue", linewidth = 3)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  ggtitle("Y basis procrustes to X basis -> curve")

sqrt(mean((as.matrix(curveX[,-1])-as.matrix(curveY[,-1]))^2))
# 0.7469886

#\\
#\\
# Align Y BASIS TO X WITH PROCRUSTES THEN CALCULATE TRAJECTORIES AND ROTATE THESE AGAIN
#\\
#\\

curveY_rot_to_curveX_procrust <- procrustes(curveX[, -1], curveY[,-1], scale = FALSE)
curveY_rot_to_curveX <- curveY_rot_to_curveX_procrust$Yrot
colnames(curveY_rot_to_curveX) <- paste0("dim",1:ncol(curveY_rot_to_curveX))

Y_basis_rotated_to_X_basis_curve_rotated_to_X_curve_p <- ggplot()+
  geom_point(data = X_reduced, aes(dim1, dim2), size = 3, color = "salmon")+
  geom_path(data = curveX, aes(dim1,dim2), color = "salmon", linewidth = 3)+
  geom_point(data = reducedY_by_rotated_v, aes(dim1, dim2), size = 3, color = "skyblue")+
  geom_path(data = curveY_rot_to_curveX , aes(dim1,dim2), color = "skyblue", linewidth = 3)+
  theme_bw()+
  theme(aspect.ratio = 1)+
  ggtitle("Y basis procrustes to X basis -> curve > Y curve procrustes to X curve")


sqrt(mean((as.matrix(curveX[,-1])-as.matrix(curveY_rot_to_curveX))^2))
# 0.3384143

all <- Y_scld_to_X_times_X_basis_p+Y_basis_rotated_to_X_basis_p+Y_basis_rotated_to_X_basis_curve_rotated_to_X_curve_p

ggsave("results/y-to-x-alignment-comparison.pdf", all, height = 10, width = 30)
