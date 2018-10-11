source('R/kpl.R')                       # kernel players
source('R/hlp.R')
source('R/kpl.R')
source("R/mnq.R")
source("R/utl.R")
source("sim/utl.R")
source("sim/mdl.R")

pcs.knl <- function(knl, psd=FALSE)
{
    N <- nrow(knl[[1]])
    r <- matrix(0, N, N)
    i <- upper.tri(r)

    ## pca
    knl <- do.call(data.frame, lapply(knl, `[`, i))
    knl <- scale(knl)
    pcs <- with(svd(knl), u %*% diag(d))

    ## reconstruct kernel matrices
    pcs <- lapply(pcs, function(f)
    {
        r[i] <- f; r <- t(r); r[i] <- f; r
    })

    ## try to make d kernels PSD
    pcs <- lapply(pcs, function(x)
    {
        v <- eigen(x, TRUE, TRUE)$values
        if(sum(v < 0) > sum(v >= 0))
            x <- -x
        x
    })

    names(pcs) <- sprintf("p%02d", seq_along(knl))
    pcs
}

knl.hom <- function(x, y, d=FALSE) cor(x[upper.tri(x, d)], y[upper.tri(y, d)])

up1 <- function(.) .[upper.tri(., TRUE)]
up2 <- function(.) .[upper.tri(., FALSE)]

dfr <- function(k)
{
    as.data.frame(do.call(cbind, lapply(k, up2)))
}

cr <- function(...)
{
    dt <- lapply(list(...), as.vector)
    dt <- do.call(cbind, dt)
    cor(dt)
}
## r <- main(N=400, P=10000, frq=0.1, mdl=~GS3+LN3)
main <- function(N=1000, P=10000, frq=.01, mdl=~AD+DM)
{
    rds <- get.rds('sim/dat')
    gls <- list(readRDS(rds))
    dat <- get.gmx(gls, N, P, Q=1, R=1)

    gmx <- with(dat, gmx[dvp > 0, ])
    knl.gmx <- krn(gmx, mdl)

    P <- ncol(gmx)
    fmx <- gmx[, sample(c(rep(TRUE, P * frq), rep(FALSE, P - P * frq)))]
    knl.fmx <- krn(fmx, mdl)
    
    ## correlation between the whole SNP and function SNP kernel
    hom <- mapply(knl.hom, knl.gmx, knl.fmx, MoreArgs=list(d=TRUE))

    dfr.gmx <- dfr(knl.gmx)
    dfr.fmx <- dfr(knl.fmx)

    cor.gmx <- round(cor(dfr.fmx), 2)
    cor.fmx <- round(cor(dfr.fmx), 2)
    list(dfr.gmx=dfr.gmx, dfr.fmx=dfr.fmx,
         knl.gmx=knl.gmx, knl.fmx=knl.fmx, hom=hom,
         cor.gmx=cor.gmx, cor.fmx=cor.fmx)
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
