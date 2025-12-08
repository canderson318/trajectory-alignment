
rm(list = ls()); gc(verbose = FALSE)

pacman::p_load(dplyr, ggplot2, tidyr, magrittr)

# nxp data
set.seed(190230)
X <- replicate(10, rnorm(100, 10, 2), simplify = "array")

# nxn identity matrix
I <- rep(0, nrow(X))%*% t(rep(0, nrow(X))) %>% 
  {diag(.) <- 1; .}

# nxn projection matrix: projects any vector onto subspace spanned by means 
P <- (rep(1, nrow(X))%*%t(rep(1, nrow(X)))) / nrow(X)

# #+ PX gives column means repolicated across all rows
# all(round((P%*%X)[1,, drop = T] , 3) == round(colMeans(X), 3))

# subtract mean from each column of X
X_cent <- (I-P) %*% X

frobenius_norm <-sqrt(sum( diag(t(X_cent)%*%X_cent) ))

X_scl <- X_cent/frobenius_norm

summary(X_scl)

