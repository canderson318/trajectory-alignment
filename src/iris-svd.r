
rm(list = ls()); gc()
pacman::p_load(dplyr, magrittr, zeallot, ggplot2, myPackage, install = FALSE)

wd = readLines("~/.wd.txt")

tmp_plot  <-  function( path = sprintf("%s/results/temp-plot.png", wd), ...){
  temp_plot(path = path, ...)
}


dat = iris
names(dat) = snakecase::to_snake_case(names(dat))

datm = data.matrix(dat[,-ncol(dat)])
datm = scale(datm)

c( s, U, V) %<-% svd(datm)
Vt = t(V)
S = diag(s)

pcs = data.frame(prcomp(datm, cale = FALSE)$x)
names(pcs) <- str_replace(names(pcs), 'PC', "X")

(
    # data.frame( U%*%S%*%Vt ) %>% # full reconstruction
    data.frame( datm%*%V ) %>% # Feature/PC Space
    # data.frame( U%*%S ) %>% # Feature/PC Space
    # data.frame( U%*%Vt) %>% # 
    # data.frame( pcs) %>% # PCs
        ggplot( aes(X1, X2))+
            geom_point(aes(color = dat$species), size = 2)+
            theme_bw()
) %>% tmp_plot(plot = .)

(
    ggplot(dat, aes(sepal_length,sepal_width))+
        geom_point(aes(color = dat$species), size = 2)+
            theme_bw()
) %>% tmp_plot(plot = .)


