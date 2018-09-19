source('R/kpl.R')                       # kernel players

source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")
source("R/utl.R")
source("sim/utl.R")
kpl.t1 <- function(N=100, P=300)
{
    x <- readRDS('data/1kg_c05.rds')[1:N, 1:P]

    k1 <- kin(x)
    k2 <- gau(x)
    k3 <- esn(x, p=1)
    k4 <- esn(x, p=.5)
    
    par(mfrow=c(2, 2))
    image(k1)
    image(k2)
    image(k3)
    image(k4)
}

knl.lm <- function(knl, coef=NULL)
{
    N <- nrow(knl[[1]])
    knl <- c(list(diag(N)), knl)

    udx <- upper.tri(knl[[1]], TRUE)
    kdx <- seq(3, l=(length(knl) - 2))

    dat <- do.call(cbind, lapply(knl, `[`, udx))

    print(cor(dat))
    for(i in kdx)
    {
        m <- lm(dat[, i] ~ dat[, 1:(i-1)])
        dat[, i] <- m$residuals
    }

    print(cor(dat))
    for(i in kdx)
    {
        k <- knl[[i]]
        k[udx] <- dat[, i]
        k <- t(k)
        k[udx] <- dat[, i]
        knl[[i]] <- k
    }
    
    for(i in kdx)
    {
        knl[[i]] <- std(knl[[i]])
    }
    knl[-1]
}

knl.pc <- function(knl)
{
    udx <- upper.tri(knl[[1]], TRUE)
    dat <- do.call(cbind, lapply(knl, `[`, udx))
    pcs <- with(svd(dat), u %*% diag(d))

    for(i in seq_along(knl))
    {
        k <- knl[[i]]
        k[udx] <- pcs[, i]
        k <- t(k)
        k[udx] <- pcs[, i]
        knl[[i]] <- k
    }
    knl
}

