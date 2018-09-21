source('R/kpl.R')                       # kernel players

source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")
source("R/utl.R")
source("sim/utl.R")

as.dfr <- function(knl)
{
    if(!is.data.frame(knl))
    {
        i <- upper.tri(knl[[1]], TRUE)
        knl <- do.call(data.frame, lapply(knl, `[`, i))
    }
    knl
}

as.kmx <- function(dfr)
{
    n <- nrow(dfr)
    n <- (-1 + sqrt(1 + 8 * n))/2
    r <- matrix(0, n, n)
    i <- upper.tri(r, TRUE)
    lapply(dfr, function(f)
    {
        r[i] <- f
        r <- t(r)
        r[i] <- f
        r
    })
}

as.mdl <- function(knl, mdl=NULL)
{
    if(is.null(mdl))
        mdl <- paste('~', paste(names(knl), collapse='+'))
    if(is.character(mdl))
        mdl <- as.formula(mdl)
    mdl
}

evl.knl <- function(knl, mdl=NULL)
{
    knl <- as.dfr(knl)
    mdl <- as.mdl(knl)

    ## terms and varialbes
    tmr <- terms(mdl)
    attr(tmr, 'intercept') <- 0
    attr(tmr, 'response') <- 0

    mmx <- as.data.frame(model.matrix(tmr, knl))
    mmx
}


pcs.knl <- function(knl)
{
    pcs <- data.frame(with(svd(as.dfr(knl)), u %*% diag(d)))
    names(pcs) <- sprintf("p%02d", seq_along(knl))

    pcs <- as.kmx(pcs)
    pcs <- lapply(pcs, function(x)
    {
        v <- eigen(x, TRUE, TRUE)$values
        if(sum(v < 0) > sum(v >= 0)) x <- -x
        x
    })
    pcs <- as.dfr(pcs)
    pcs
}


main <- function(N=1000, P=10000, frq=.01, oks=c(ad, dm, rs, ht, ga, p1))
{
    rds <- get.rds('sim/dat')
    gls <- list(readRDS(rds))
    dat <- get.gmx(gls, N, P, Q=1, R=1)

    gmx <- with(dat, gmx[dvp > 0, ])
    ## gmx <- scale(gmx)

    knl.gmx <- krn(gmx, oks)
    knl.gmx <- evl.knl(knl.gmx, mdl=NULL)

    P <- ncol(gmx)
    fmk <- sample(c(rep(TRUE, P * frq), rep(FALSE, P - P * frq)))
    fmx <- gmx[, fmk]
    knl.fmx <- krn(fmx, oks)
    knl.fmx <- evl.knl(knl.fmx, mdl=NULL)
    
    ## correlation between the whole SNP and function SNP kernel
    cor.fgx <- sapply(colnames(knl.gmx), function(n) cor(knl.fmx[, n], knl.gmx[, n]))

    ## pca
    pcs.gmx <- pcs.knl(knl.gmx)
    pcs.fmx <- pcs.knl(knl.fmx)
    
    list(knl.gmx=knl.gmx, knl.fmx=knl.fmx,
         cor.fgx=cor.fgx,
         pcs.gmx=pcs.gmx, pcs.fmx=knl.fmx)
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

