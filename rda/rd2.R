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

mix.knl <- function(kns, std=1, ...)
{
    k <- Reduce(`*`, kns)
    if(is.na(std))
        std <- 1
    if(std)
        k <- sign(k) * abs(k)^(1 / length(kns))
    k <- list(k)
    names(k) <- paste0(names(kns), collapse=":")
    k
}

get.knl <- function(pfx, bat, rnd=~uwt, std=1, ...)
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
        mix.knl(kns[fc[, f] > 0], std)[[1]]
    },
    simplify=FALSE)
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
main <- function(mdl=sbp~sex+age, rnd=~sgn, ply=1, pfx='rda/gwk', bat='030', mtd=fwd, ...)
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
    
    ## sample matching
    . <- Reduce(intersect, lapply(c(list(y, X), K), rownames))
    y <- y[., , drop=FALSE]
    X <- X[., , drop=FALSE]
    K <- lapply(K, `[`, ., .)
    
    ## level 1 model development
    ft1 <- mtd(y, K, X, ...)
    vc1 <- intersect(names(K), names(ft1$par))

    ## level 2
    vc2 <- setdiff(names(K), vc1)
    if(ply > 1 && length(vc2) > 0)
    {
        KN1 <- K[vc1]
        KN2 <- K[vc2]
        K <- KN1
        for(k1 in names(KN1))
        {
            KNX <- list()
            for(k2 in names(KN2))
            {
                KNX <- c(KNX, mix.knl(c(KN1[k1], KN2[k2]), ...))
            }
            md2 <- mtd(y, KNX, X, mtd=fwd, msl=1)
            vc2 <- intersect(names(KNX), names(md2$par))
            K <- c(K, KNX[vc2])
        }
        ft2 <- mtd(y, K, X, ...)

        if(ft2$rpt$nlk < ft1$rpt$nlk &&
           ft2$rpt$yel < ft1$rpt$yel &&
           ft2$rpt$ycl > ft1$rpt$ycl)
            fit <- ft2
        else
            fit <- ft1
    }
    else
        fit <- ft1
    ## return
    r <- within(fit, arg <- arg)
    r
}

colMedian <- function(x) apply(x, 2L, median)

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
    
    ## perform testing for each model configuration
    ## evs <- lapply(split(seq.int(nrow(ags)), ags), function(gdx)
    evs <- by(pas, ags, function(grp)
    {
        ## unique arguments and aggregated configuration
        arg <- unique(grp[, names(ags)])
        par <- grp[grp$bat %in% mdx, grep("X00", names(grp), TRUE):ncol(grp)]
        par[is.na(par)] <- 0.0
        par <- switch(agg, md=colMedian, xt=agg.vcs, colMeans)(par)
        fxs <- fx(par)
        vcs <- vc(par)
        vcs <- vcs[vcs > 0]
        rnd <- paste("~", paste(names(vcs)[-1], collapse="+"))
        arg$rnd <- gsub('[.]', ':', rnd)
        arg$bat <- edx
        arg$nbt <- nbt
        arg$agg <- agg
        arg$evl <- "avg"
        
        ## response variable and fix effect design matrix 
        D <- get.dsg(arg$mdl, phe, pcs)
        
        ## for evaluation: fetch batch kernel, polynomial expansion, standardize
        K <- do.call(get.knl, arg)      # NOTE: arg$bat == edx !!
        
        ## matching
        i <- Reduce(intersect, lapply(c(K, D), rownames))
        y <- D$y[i, , drop=FALSE]
        X <- cbind(D$X[i, ], D$pcs[i, ])
        K <- lapply(K, `[`, i, i)
        
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


asmb <- function(mds, nbt=Inf, edx=NULL, sbt=NULL, agg="mu", nes=Inf)
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

    ## number of ensembles
    nes <- max(1, min(nes, length(bdx)))
    
    ## perform testing for each model configuration
    ## evs <- lapply(split(seq.int(nrow(ags)), ags), function(gdx)
    evs <- by(pas, ags, function(grp)
    {
        ## unique arguments and aggregated configuration
        arg <- unique(grp[, names(ags)])
        par <- grp[grp$bat %in% mdx, grep("X00", names(grp), TRUE):ncol(grp)]
        par[is.na(par)] <- 0.0
        fag <- switch(agg, md=colMedian, xt=agg.vcs, colMeans)
        mpa <- fag(par)
        vcs <- vc(mpa)
        vcs <- vcs[vcs > 0]
        rnd <- paste("~", paste(names(vcs)[-1], collapse="+"))
        arg$rnd <- gsub('[.]', ':', rnd)
        arg$bat <- edx
        arg$nbt <- nbt
        arg$agg <- agg
        arg$evl <- "esb"
        
        ## response variable and fix effect design matrix 
        D <- get.dsg(arg$mdl, phe, pcs)
        ## for evaluation: fetch batch kernel, polynomial expansion, standardize
        K <- do.call(get.knl, arg)      # NOTE: arg$bat == edx !!
        
        ## matching
        i <- Reduce(intersect, lapply(c(K, D), rownames))
        y <- D$y[i, , drop=FALSE]
        X <- cbind(D$X[i, ], D$pcs[i, ])
        K <- lapply(K, `[`, i, i)
        
        ## rpt <- rbind(rpt, hat)
        hat <- lapply(1:nes, function(i)
        {
            print(i)
            j <- sample.int(nrow(par), nrow(par) * 0.5, replace=TRUE)
            w <- fag(par[j, ])
            r <- loo(y, K, X, w)
            r
        })
        hat <- do.call(cbind, hat)
        hat <- rowMeans(hat)

        yel <- mean(y - hat)
        ycl <- cor(y, hat)
        ret <- list(arg=arg, rpt=data.frame(yel=yel, ycl=ycl), par=mpa)
        ret
    })
    
    ## pack, round up, and return
    arg <- do.call(rbind, evs %$% 'arg')
    par <- do.call(rbind, evs %$% 'par')
    rpt <- do.call(rbind, evs %$% 'rpt')
    
    ret <- list(arg=arg, par=par, rpt=rpt)
    ret
}
