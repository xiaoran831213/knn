library(MASS)
library(plinkBED)
library(wrapGCTA)
source('R/hlp.R')                       # helpers
source('R/kpl.R')                       # kernel players
source('R/utl.R')                       # utilities
source('R/vcm.R')                       # variance component models (VCM)
source('R/mnq.R')                       # MINQUE
source('R/ply.R')
source('rda/rpt.R')
source('rda/abc.R')

scn <- function(., w) scan(., w, sep='\t', skip=1, quiet=TRUE, blank.lines.skip=FALSE)

get.phe <- function(use.cache=TRUE)
{
    rds <- 'rda/phe.rds'
    if(file.exists(rds) && use.cache)
        phe <- readRDS(rds)
    else
    {
        ## sample ID
        eid <- scn('rda/phe/eid.txt', 0L)
        nsb <- length(eid)

        ## covariate
        cvr <- within(list(),
        {
            sex <- scn('rda/phe/sex.txt', 0L)
            age <- scn('rda/phe/age.txt', 0L)
            age <- scale(age)

            ## smk <- scn('rda/phe/smk_bin.txt', 0L)
            sbp <- scn('rda/phe/sbp.txt', 0L)
            sbp <- log(1 + sbp)
            
            dbp <- scn('rda/phe/dbp.txt', 0L)
            dbp <- log(1 + dbp)

            smk <- scn('rda/phe/smk_frq.txt', 0.0)
            alc <- scn('rda/phe/alc_frq.txt', 0.0)
            cnb <- scn('rda/phe/cnb_frq.txt', 0.0)

            smk <- scale(smk)
            alc <- scale(alc)
            cnb <- scale(cnb)
        })
        vol <- get.vol(FALSE)

        ## combine and return
        phe <- do.call(data.frame, c(cvr, vol))
        rownames(phe) <- eid
        saveRDS(phe, rds)
    }
    invisible(phe)
}

## histogram
get.hist <- function(out="~/img/rda_mix_hst.png", use.cache=TRUE)
{
    options(stringsAsFactors=FALSE)

    rds <- 'rda/hst.rds'
    if(file.exists(rds) && use.cache)
        dat <- readRDS(rds)
    else
    {
        eid <- scn('rda/phe/eid.txt', 0L)
        dat <- within(list(),
        {
            sex <- scn('rda/phe/sex.txt', 0L)
            age <- scn('rda/phe/age.txt', 0L)

            smk <- scn('rda/phe/smk_frq.txt', 0.0)
            alc <- scn('rda/phe/alc_frq.txt', 0.0)
            cnb <- scn('rda/phe/cnb_frq.txt', 0.0)
        })
        dat <- do.call(data.frame, dat)
        rownames(dat) <- eid

        dat <- melt(dat, variable.name='rsp', value.name='val')
        dat <- subset(dat, !is.na(val))
        dat <- rbind(cbind(dat, zero='with 0'), cbind(subset(dat, val!=0), zero='without'))
        saveRDS(dat, rds)
    }

    d <- subset(dat, rsp %in% c('smk', 'cnb', 'alc'))
    d <- within(d,
    {
        val <- pmin(val, quantile(val, .995))
        ## val  <- pmax(val, quantile(val, .025))
        grp <- paste(rsp, zero)
    })
    
    g <- ggplot(d)
    ## r <- 6
    ## g <- g + xlim(-r, +r) + ylim(-r, +r)
    ## g <- g + ylab(quote(bold(eta) %~% N(0, Sigma)))
    ## g <- g + xlab(quote(tilde(bold(eta))))

    ## histogram plot in the back in gray
    ## g <- g + geom_histogram(aes(x=val, y=..density.., fill=rsp), d, alpha=1, binwidth=.5)
    g <- g + geom_histogram(aes(x=val, fill=rsp), d, color='black', bins=10, binwidth=.5)
    g <- g + facet_wrap(~grp, nrow=2, scales='free', dir='v')
    g <- g + theme(
        legend.position = "none",
        strip.text=element_text(face='bold', size=15),
        strip.background = element_rect(colour="red", fill="#CCCCFF"),
        axis.title = element_text(face='bold', size=15))

    ggsave(out, width=19, height=7)
    invisible(dat)
}

get.pcs <- function(use.cache=TRUE)
{
    pcs <- 'rda/phe/gno_pcs.txt'
    eid <- 'rda/phe/eid.txt'
    rds <- sub('[.].*$', '.rds', pcs)
    if(file.exists(rds) && use.cache)
        pcs <- readRDS(rds)
    else
    {
        eid <- scn(eid, 0L)
        pcs <- scn(pcs, .0)
        pcs <- matrix(pcs, length(eid), byrow=TRUE)
        colnames(pcs) <- sprintf('p%02d', seq.int(ncol(pcs)))
        rownames(pcs) <- eid
        pcs <- as.data.frame(pcs)
        saveRDS(pcs, rds)
    }
    invisible(pcs)
}

get.vol <- function(use.cache=TRUE)
{
    vol <- 'rda/phe/vol.txt'
    eid <- 'rda/phe/eid.txt'
    rds <- sub('[.].*$', '.rds', vol)
    if(file.exists(rds) && use.cache)
        vol <- readRDS(rds)
    else
    {
        eid <- scn(eid, 0L)
        vol <- scn(vol, .0)
        vol <- matrix(vol, length(eid), byrow=TRUE)
        colnames(vol) <- sprintf('v%02d', seq.int(0L, ncol(vol) - 1L))
        rownames(vol) <- eid
        vol <- as.data.frame(scale(vol))
        saveRDS(vol, rds)
    }
    invisible(vol)
}


get.knl <- function(pfx, bat, rnd=~cov, ply=2, std=0, ...)
{
    kns <- sapply(all.vars(as.formula(rnd)), function(k)
    {
        if(k == 'rel' || k == 'std')
            ret <- readREL(paste0(file.path(pfx, bat), paste0('.', k)))
        if(k == 'ibs')
            ret <- readIBS(paste0(file.path(pfx, bat), paste0('.', k)))
        nas <- sort(unique(as.vector(which(is.na(ret), TRUE))))
        if(length(nas) > 0L)
            ret <- ret[-nas, -nas]
        ret
    }, simplify=FALSE)
    kns <- PK(kns, D=ply)()
    
    if(std)
        kns <- lapply(kns, function(k) k / mean(diag(k)))
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

get.dsg <- function(mdl, phe, pcs)
{
    gus <- rnorm(nrow(phe))
    psi <- rnorm(nrow(phe), 1) - 1
    bin <- rbinom(nrow(phe), 4, .5)
    ## raw data
    dat <- cbind(phe, pcs, gus=gus, psi=psi, bin=bin)
    
    ## model matrix
    mfr <- model.frame(mdl, dat)

    ## response variable
    y <- as.matrix(model.response(mfr))
    colnames(y) <- names(mfr)[1]

    ## design matrix: conventional covariates
    X <- model.matrix(mfr, dat)[, -1]
    ## colnames(X)[1] <- 'X00'
    
    ## principle components
    ipc <- names(mfr)[names(mfr) %in% names(pcs)]
    if(length(ipc) > 0)
    {
        i <- intersect(rownames(X), rownames(pcs))
        X <- cbind(X[i, !colnames(X) %in% ipc, drop=FALSE],
                   as.matrix(pcs[i, names(pcs) <= max(ipc), drop=FALSE]))
    }

    list(y=y, X=X)
}

MQ0 <- function(y, K, X, itr=50, ...) knl.mnq(y, K, X, zbd=0, vbd=0, ebd=0, itr=itr, ...)
MQ1 <- function(y, K, X, itr=50, ...) knl.mnq(y, K, X, zbd=0, vbd=1, ebd=0, itr=itr, ...)
MQ2 <- function(y, K, X, itr=50, ...) knl.mnq(y, K, X, zbd=0, vbd=0, ebd=1, itr=itr, ...)
MQ3 <- function(y, K, X, itr=50, ...) knl.mnq(y, K, X, zbd=0, vbd=1, ebd=1, itr=itr, ...)

ZQ0 <- function(y, K, X, itr=50, ...) knl.mnq(y, K, X, zbd=1, vbd=0, ebd=0, itr=itr, ...)
ZQ1 <- function(y, K, X, itr=50, ...) knl.mnq(y, K, X, zbd=1, vbd=1, ebd=0, itr=itr, ...)
ZQ2 <- function(y, K, X, itr=50, ...) knl.mnq(y, K, X, zbd=1, vbd=0, ebd=1, itr=itr, ...)
ZQ3 <- function(y, K, X, itr=50, ...) knl.mnq(y, K, X, zbd=1, vbd=1, ebd=1, itr=itr, ...)

ZGC <- function(y, K, X, itr=100, ...) gcta.reml(y, K, X[, -1], itr=itr, zbd=1, alg=1)
UGC <- function(y, K, X, itr=100, ...) gcta.reml(y, K, X[, -1], itr=itr, zbd=0, alg=1)

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
main <- function(mdl=alc~sex+age+p02, rnd=~std, ply=2, pfx='rda/wkn/201', bat='030', mtd=MQ3, ...)
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
    ## K <- get.knl(pfx, bat, rnd, ply, std=1, ...)
    K <- cal.knl(pfx, bat, rnd, ply, ...)
    
    ## sample matching
    . <- Reduce(intersect, lapply(c(list(y, X), K), rownames))
    y <- y[., , drop=FALSE]
    X <- X[., , drop=FALSE]
    K <- lapply(K, `[`, ., .)
    
    ## model development
    r <- mtd(y, K, X, ...)
    r <- within(r, arg <- arg)
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
        par <- switch(agg, md=colMedian, xt=agg.vcs, colMeans)(par)
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

        rpt <- vpd(na.omit(par), y, K, X, rt=0)
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

ts1 <- function(v=c(1.0, 1.0), s=v, pmu=0)
{
    N <-  500
    P <- 5000
    ## K <- get.knl('rda/wkn/201', '006', ~std, ply=1)[[1]]
    X <- scale(matrix(rnorm(N * P, sd=2), N, P))
    K <- tcrossprod(X) / P
    
    ## generating
    N <- nrow(K)
    y <- drop(MASS::mvrnorm(1, rep(0, N), v[2] * K)) + rnorm(N, sd=sqrt(v[1]))

    ## evaluate
    if(pmu)
        y <- sample(y)

    F <- s[2] * K
    E <- diag(s[1], N)
    V <- E + F
    A <- solve(V)
    Ay <- A %*% y
    
    blp <- drop(F %*% Ay)
    loo <- drop(y - Ay / diag(A))

    kpa <- c(F=kappa(F), V=kappa(V), FV=kappa(F %*% A))
    print(list(kpa=kpa))

    slp <- c(coef(lm(y ~ blp))[2], coef(lm(y ~ loo))[2])
    print(list(SLP=slp))
    
    crr <- c(blp=cor(y, blp), loo=cor(y, loo))
    print(list(COR=crr))
}

## program Jian Yang's GRM kernel.
## make sure swt() is equivalent to GCTA GRM.
ts2 <- function()
{
    gmx <- readBED("rda/vkn/000")
    std <- unname(readREL("rda/vkn/000")[1:10, 1:10])

    q <- colMeans(gmx, na.rm=TRUE) / 2  # AF
    k <- swt(gmx[1:10, ], o=1, q=q)

    print(all.equal(k, std))
}


ts3 <- function(r=100, gmx=NULL, out=NULL)
{
    if(is.character(gmx))
    {
        gmx <- readRDS(gmx)
    }
    q <- colMeans(gmx, na.rm=TRUE) / 2  # AF
    N <- nrow(gmx)
    P <- ncol(gmx)
    ret <- list()
    for(r in 1:r)
    {
        cat("run=", r, "\r", sep="")
        x <- gmx[sample.int(N, 2), ]
        for(p in c(2e2, 2e3, 2e4, 2e5))
        {
            j <- sample.int(P, p)
            k <- swt(x[, j], o=1, q=q[j])
            ret[[length(ret) + 1]] <- DF(idx=r, P=p, K=k[1, 2])
        }
    }
    cat("\n")
    ret <- do.call(rbind, ret)
    ret <- within(ret, P <- factor(P))

    g <- ggplot(rp3, aes(x=P, y=K))
    g <- g + geom_boxplot()
    g <- g + coord_cartesian(ylim=c(-.2, +.2))

    if(!is.null(out))
        ggsave(out, g)
    g
}

ts4 <- function(..., out=NULL)
{
    fs <- list(...)
    d0 <- lapply(fs, function(f)
    {
        r <- readREL(f)
        r <- r[upper.tri(r)]
        data.frame(x=sub("^.*[/]", "", f), y=r)
    })
    d0 <- do.call(rbind, d0)

    g <- ggplot(d0, aes(x=x, y=y))
    g <- g + geom_boxplot()
    g <- g + coord_cartesian(ylim=c(-.2, +.2))

    if(!is.null(out))
        ggsave(out, g)
    g
}
