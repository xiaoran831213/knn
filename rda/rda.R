library(MASS)
library(plinkBED)
library(wrapGCTA)
library(mnq)
source('R/hlp.R')                       # helpers
source('R/kpl.R')                       # kernel players
source('R/utl.R')                       # utilities
source('R/vcm.R')                       # variance component models (VCM)
source('R/ply.R')
source('rda/rpt.R')
source('rda/abc.R')
source('rda/phe.R')

get.knl <- function(pfx, bat, rnd=~uwt, ply=2, std=0, ...)
{
    fns <- dir(pfx, bat, full=TRUE)
    dat <- sub(paste0("^.*", bat, "[.]"), "", fns)
    names(dat) <- dat
    dat <- data.frame(t(dat), stringsAsFactors=FALSE)[0, ]
    tm <- attributes(terms(as.formula(rnd), data=dat))
    fc <- tm$factors[rowSums(tm$factors) > 0, , drop=FALSE]

    ## read basic terms
    kns <- sapply(rownames(fc), function(f)
    {
        print(paste0(file.path(pfx, bat), ".", f))
        k <- readRDS(paste0(file.path(pfx, bat), ".", f))
        n <- sort(unique(as.vector(which(is.na(k), TRUE))))
        if(length(n > 0))
            k <- k[-n, -n]
        k
    }, simplify=FALSE)
    iid <- Reduce(intersect, lapply(kns, rownames))
    kns <- lapply(kns, `[`, iid, iid)
    
    ## expand terms
    kns <- sapply(colnames(fc), function(f)
    {
        k <- Reduce(`*`, kns[fc[, f] > 0])
    }, simplify=FALSE)
    kns
}

cal.knl <- function(pfx, bat, rnd=~S1T, ply=1, ...)
{
    x <- readBED(paste0(file.path(pfx, bat)))
    i <- readFAM(paste0(file.path(pfx, bat)))
    rownames(x) <- i$iid
    
    kns <- sapply(all.vars(as.formula(rnd)), function(k)
    {
        do.call(k, list(x))
    },
    simplify=FALSE)

    kns <- PK(kns, D=ply)()
    kns
}

dvs <- function(K, r=100)
{
    if(is.matrix(K))
    {
        K <- list(K)
        single <- TRUE
    }
    else
        single <- FALSE

    n <- nrow(K[[1]])
    z <- matrix(rnorm(n * r), n, r)
    v <- lapply(K, function(k)
    {
        mean(z * crossprod(k, z))
    })

    if(single)
        v <- v[[1]]
    v
}
    
#' analyze a real data
main <- function(mdl=sbp~sex+age, rnd=~sgn, ply=1, pfx='rda/gwk', bat='030', mtd=fcv, ...)
{
    options(stringsAsFactors=FALSE)
    arg <- as.list(match.call(expand.dots=TRUE)[-1])
    idx <- !sapply(arg, is.vector)
    arg[idx] <- lapply(arg[idx], deparse)
    arg <- do.call(data.frame, arg)
    if(!is.null(arg$seed) && arg$seed != "NULL")
        set.seed(arg$seed)
    pmu <- arg$pmu %||% 0

    ## make design matrices for fixed effect and response
    . <- get.dsg(mdl, get.phe(), get.pcs())
    X <- .$X
    y <- .$y
    if(pmu)
        rownames(y) <- sample(rownames(y))
    ## read kernels, expand, and maybe standardize
    K <- get.knl(pfx, bat, rnd, ply, ...)
    ## K <- cal.knl(pfx, bat, rnd, ply, ...)
    
    ## sample matching
    . <- Reduce(intersect, lapply(c(list(y, X), K), rownames))
    y <- y[., , drop=FALSE]
    X <- X[., , drop=FALSE]
    K <- lapply(K, `[`, ., .)
    
    ## level 1 model development
    md1 <- mtd(y, K, X, ...)
    
    ## prepare level 2:
    r <- within(md1, arg <- arg)
    r
}

#' aggregate the estimated variance component
agg.vcs <- function(pas)
{
    cfg <- pas[, c('mdl', 'mtd')]
    pas <- pas[, grep('X00', names(pas), TRUE):ncol(pas)]
    agg <- by(data=pas, INDICES=cfg, FUN=function(par)
    {
        nes <- nrow(par)                    # number of estimates
        npa <- ncol(par)
        rownames(par) <- sprintf("%04x", seq.int(nes))

        idx <- grep("EPS", names(par), TRUE)
        vcs <- seq(idx, npa)                # variance components

        eps <- par[, idx]                   # noise
        emd <- median(eps)                  # median of noise
        
        ## upper half and lower half according to noise
        uvc <- par[eps > emd, ]
        lvc <- par[eps < emd, ]
        uvc <- uvc[order(uvc[, idx], decreasing=FALSE), ]
        lvc <- lvc[order(lvc[, idx], decreasing=TRUE ), ]
        nes <- min(nrow(uvc), nrow(lvc))
        uvc <- uvc[seq(nes), ]
        lvc <- lvc[seq(nes), ]
        ## the median, can be empty when nes is even, or multiple lines.
        mvc <- par[eps == emd, , drop=FALSE]
        
        ## bi-directional cumulative mean, starting from the median
        cmb <- (uvc + lvc) / 2
        if(nrow(mvc) > 0)
            cmb <- rbind(cmb, na.omit(colMeans(mvc)))

        cmm <- do.call(data.frame, lapply(cmb, cummean)) ## [-1, ]
        
        ## number of negatives
        neg <- rowSums(cmm[, vcs] <  0, na.rm=TRUE)

        ## level of sparsity
        spa <- rowSums(cmm[, vcs] == 0, na.rm=TRUE)
        
        ## most sparse, non-negative, and inclusive solution
        msk <- neg == min(neg)
        ## msk <- spa == max(spa)
        ret <- cmm[msk, ]
        ret <- unlist(ret[nrow(ret), , drop=TRUE])
        ret[idx] <- emd
        ret
    })
    agg <- cbind(expand.grid(dimnames(agg)), do.call(rbind, agg))
    agg
}


## p1 <- plt1('rda/dvp/201_gus')
## d2 <- subset(p1, mdl=='gus~(std+ibs)^2' & mtd=='ZMQ' & key=='EPS' & val < .1)
#' cross valication
#'
#' @param mds a directory of models
#' @param ncv integer number cross validation.
#' @param edx evaluation data index
cval <- function(mds, nbt=Inf, edx=NULL, sbt=NULL, agg="mu")
{
    cfg <- as.list(match.call(expand.dots=TRUE)[-1])
    sbt <- substitute(sbt)
    
    ## phenotypes
    phe <- get.phe()
    pcs <- get.pcs()
    
    ## model paramters, configurations, and reports
    ags <- agg.arg(mds, simple=FALSE)
    pas <- agg.par(mds, with.arg=TRUE, simple=FALSE)
    rps <- agg.rpt(mds, with.arg=TRUE, simple=FALSE)
    
    ## select models by sbt
    if(is.call(sbt))
    {
        msk <- eval(sbt, ags, parent.frame())
        if(length(msk) == 0)
            return(NULL)
        pas <- pas[msk, ]
        rps <- rps[msk, ]
        ags <- ags[msk, ]
    }
    
    ## handle batches
    bts <- ags$bat; ags$bat <- NULL

    ## find batches universally available for all configurations
    ags[is.na(ags)] <- 0
    bdx <- Reduce(intersect, split(bts, ags))

    ## select one batch for testing, and {nbt} batches for modeling
    if(is.null(edx) || !edx %in% bdx)   # one for testing
        edx <- sample(bdx, 1L)
    ## from the rest, choose {nbt} to create an average model
    nbt <- max(1, min(nbt, length(bdx) - 1))
    mdx <- sort(sample(setdiff(bdx, edx), nbt))
    
    colMedian <- function(x) apply(x, 2L, median)
    ## perform testing for each model configuration
    ## evs <- lapply(split(seq.int(nrow(ags)), ags), function(gdx)
    evs <- by(pas, ags, function(grp)
    {
        ## unique arguments and aggregated configuration
        arg <- unique(grp[, names(ags)])
        par <- grp[grp$bat %in% mdx, grep("X00", names(grp), TRUE):ncol(grp)]
        par[is.na(par)] <- 0.0
        par <- switch(agg, md=colMedian, xt=agg.vcs, colMeans)(par)
        vcs <- par[seq(which(names(par) == "EPS"), length(par))]
        fix <- par[setdiff(names(par), names(vcs))]
        vcs <- vcs[vcs > 0]
        arg$rnd <- paste("~", paste(names(vcs)[-1], collapse="+"))
        arg$bat <- edx
        arg$nbt <- nbt
        arg$agg <- agg
        
        ## response variable and fix effect design matrix 
        D <- get.dsg(arg$mdl, phe, pcs)
        
        ## for evaluation: fetch batch kernel, polynomial expansion, standardize
        K <- do.call(get.knl, arg)      # NOTE: arg$bat == edx !!
        
        ## matching
        i <- Reduce(intersect, lapply(c(K, D), rownames))
        y <- D$y[i, , drop=FALSE]
        X <- cbind(D$X[i, ], D$pcs[i, ])
        K <- lapply(K, `[`, i, i)
        
        ## evaluate now
        ## evl <- vpd(y, K, na.omit(par), X, rt=0, hat=TRUE)
        ## rpt <- evl$rpt
        ## .x <- evl$hat$zu
        ## .y <- unname(drop(y - evl$hat$xb))
        ## hat <- kms(.x, .y, 10)
        ## hat <- within(hat, key <- as.character(key))
        ## rpt <- rbind(rpt, hat)
        rpt <- vpd(y, K, X, par)
        ret <- list(arg=arg, rpt=rpt, par=par)
        ret
    })
    
    ## pack, round up, and return
    arg <- do.call(rbind, evs %$% 'arg')
    par <- do.call(rbind, evs %$% 'par')
    rpt <- do.call(rbind, evs %$% 'rpt')
    
    ret <- list(arg=arg, par=par, rpt=rpt)
    ret
}

## diagnostic plot
pl1 <- function(dat, key, val, ...)
{
    library(ggplot2)
    mx <- key[dat == 'kmu']
    my <- val[dat == 'kmu']
    s1 <- val[dat == 'ksd']
    alpha <- 1 - 0.8 * (s1 - min(s1)) / (max(s1) - min(s1))
    cl <- rgb(0, 0, 1, alpha)
    symbols(mx, my, circles=s1, inches=1/3, ann=TRUE, bg=cl, fg=NULL, cex.main=.8, ...)
}


## r <- main(mdl=alc~sex+age, pfx="rda/wkn/201", mtd=UGC, itr=100, bat="002", ply=2, rnd=~ibs)
kms <- function(x, y, K=10)
{
    x <- unname(x)
    y <- unname(y)

    ## clusters of x, and clustering of y accordingly
    kx <- kmeans(x, K)
    xx <- kx$centers
    yx <- unname(sapply(split(y, kx$cluster), mean))
    ys <- unname(sapply(split(y, kx$cluster), sd))
    ix <- sort(xx, index.return=TRUE)$ix
    xx <- round(xx[ix], 3)              # cluster of X
    yx <- round(yx[ix], 3)              # cluster of y by x
    ys <- round(ys[ix], 3)              # cluster of y's sd
    rbind(DF(dat='kmu', key=xx, val=yx), DF(dat='ksd', key=xx, val=ys))

    ## clusters of y, and clustering of x accordingly
    ## ky <- kmeans(y, K)                  
    ## yy <- ky$centers
    ## xy <- unname(sapply(split(x, ky$cluster), mean))
    ## xs <- unname(sapply(split(x, ky$cluster), sd))
    ## iy <- sort(yy, index.return=TRUE)$ix
    ## xy <- round(xy[iy], 3)
    ## xs <- round(xs[iy], 3)
    ## yy <- round(yy[iy], 3)

    ## mm <- DF(mx=x, my=y)
    ## mm <- kmeans(mm, K)
    
    ## rbind(DF(dat='kmu', key=round(mm$centers[, 1], 3), val=round(mm$centers[, 2],3)),
    ##       DF(dat='ksd', key=round(mm$centers[, 1], 3), val=round(mm$withinss, 3)))
}


## r <- main(alc~sex+age, rnd=~c03+c08+c09, ply=1, pfx='rda/ckn', bat='025')
## r <- main(alc~sex+age, rnd=~c03+c08+c09+c10+c11, ply=1, pfx='rda/ckn', bat='025')
