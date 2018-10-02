## weight schemes

get.wgt <- function(x, m=0)
{
    p <- colMeans(x, na.rm=TRUE) / 2
    q <- 1 - p

    if(m==1)                            # beta(1, 25)
        w <- dbeta(p, 1, 25)
    else if(m==2)                       # sd (gcta)
        w <- sqrt(1 / (p * q))
    else if(m==3)                       # log
        w <- sqrt(-log10(p))
    else                                # no weight
        w <- rep(1, ncol(x))
    sweep(x, 2L, w, '*')
}

W0 <- function(x) list(W0=ply(get.wgt(x, 0))) # no weight
W1 <- function(x) list(W1=ply(get.wgt(x, 1))) # beta weight
W2 <- function(x) list(W2=ply(get.wgt(x, 2))) # sd weight
W3 <- function(x) list(W3=ply(get.wgt(x, 3))) # log weight
