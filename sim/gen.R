library(mnq)
## read genome data, simulate response variable


#' randomly pick some RDS file from one or more directories.
#'
#' @param ... the directories
#' @param n   the number of files to select
get.rds <- function(..., n=1)
{
    ds <- unlist(list(...))
    fs <- sapply(ds, dir, "[.]rds", TRUE, TRUE)
    sample(fs, n)
}

#' get and divide a segment from the genome
#' 
#' @param gmx the genomic matrix
#' @param N the size of each cohort
#' @param P the number of variants
#' @param Q the number of training cohorts
#' @param R the number of testing cohort
get.gmx <- function(gls, N=100, P=50, Q=4, R=1)
{
    ## masks
    dvp <- c(rep(1L:Q, each=N), rep(0L, N * R))
    evl <- c(rep(0L, N * Q), rep(seq(1L, l=R), each=N))
    gsp <- list()

    ## select variants and samples
    for(i in seq_along(gls))
    {
        gmx <- gls[[i]]

        ## indices
        idx <- sample.int(nrow(gmx))[seq.int(N * Q + N * R)]
        jdx <- seq(sample.int(ncol(gmx) - P, 1L), l=P)
        gmx <- gmx[idx, jdx]

        ## remove degeneracy
        af <- colMeans(gmx) / 2.0
        gmx <- gmx[, pmin(af, 1 - af) >= 0.05]
    
        ## divide
        gls[[i]] <- gmx
        gsp[[i]] <- rep(i, ncol(gmx))
    }

    ## pack up and return
    gsp <- unlist(gsp)
    gmx <- do.call(cbind, gls)
    list(gmx=gmx, dvp=dvp, evl=evl, gsp=gsp)
}

#' get variance components randomly
#'
#' @param obj the number of compoents, or a list
#' @param dfs degree of freedom for each component (rotated)
get.vcs <- function(obj, dfs=1)
{
    if(is.numeric(obj) && length(obj) == 1)
        nvc <- obj
    else
        nvc <- length(obj)
    dfs <- rep(dfs, l=nvc)
    rchisq(nvc, dfs)
}

#' get functional mask
get.fmk <- function(obj, frq=.5)
{
    if(is.numeric(obj) && length(obj) == 1)
        P <- obj
    else
        P <- NCOL(obj)
    sample(c(rep(TRUE, P * frq), rep(FALSE, P - P * frq)))
}

#' get simulated response for once
#'
#' @param gmx N x P genomic data matrix
#' @param oks formula of output kernels
#' @param eps noise size, or functional
#' @param vcs true variance components
#' @param fix true fixed coefficients
#' @param ... additional parameters
#'     dks formula of default working kernels
get.one <- function(gmx, oks=~LN, frq=.1, vcs=1, fix=0, ...)
{
    N <- nrow(gmx)
    dot <- list(...)
    lnk <- if(is.null(dot$lnk)) NL else dot$lnk
    LNK <- if(is.null(dot$LNK)) NL else dot$LNK

    ## functional variants and kernels
    fmk <- dot$fmk %||% get.fmk(gmx, frq) # functional mask
    fnl <- krn(gmx[, fmk], oks)           # functional kernel

    ## make default kernels
    if(!is.null(dot$dks))
        knl <- krn(gmx, dot$dks)
    
    ## variance components
    vcs <- rep(vcs, l=1 + length(fnl))        # propagate
    names(vcs) <- c("EPS", names(fnl))

    ## covariance matrix
    cmx <- cmb(fnl, vcs[-1], TRUE)
    ## random effect
    eta <- drop(MASS::mvrnorm(1, rep(0, N), cmx))


    ## noise function
    efn <- dot$efn %||% EGS

    ## fixed effect and mean
    names(fix) <- sprintf('X%02d', seq_along(fix) - 1)
    if(length(fix) > 1)
    {
        xmx <- replicate(length(fix) - 1, rnorm(N))
        colnames(xmx) <- names(fix[-1])
        rownames(xmx) <- rownames(gmx)
        mu <- fix[1] + xmx %*% fix[-1]
    }
    else
    {
        xmx <- NULL
        mu <- fix[1]
    }
    
    ## response = fix + rnd + eps
    rsp <- mu + LNK(lnk(eta) + efn(N, vcs[1]))
    rsp <- as.matrix(rsp)
    rownames(rsp) <- rownames(gmx)

    par <- c(fix, vcs)
    ret <- list(gmx=gmx, rsp=rsp, fmk=fmk, par=par, xmx=xmx, oks=oks)
    if(!is.null(dot$dks))
        ret <- c(ret, list(knl=knl))
    ret
}

#' write one simulated data to text
write.one <- function(ret, out)
{
    dir.create(out, FALSE, TRUE)
    rdm <- "A GENOMIC SIMULATION INSTANCE\n"

    ## genomic matrix
    gmx <- ret$gmx
    write(t(gmx), file.path(out, 'gmx.txt'), ncol(gmx), sep='\t')
    write(rownames(gmx), file.path(out, 'gmx.row'), 1)
    write(colnames(gmx), file.path(out, 'gmx.col'), 1)
    rdm <- c(rdm, "\ngmx.txt: N x P genotype matrix in dosage format.\n")
    rdm <- c(rdm, "\tgmx.row: row names of gmx, i.e., individual IDs.\n")
    rdm <- c(rdm, "\tgmx.col: column names of gmx, i.e., variant IDs.\n")
    
    ## response variable
    rsp <- ret$rsp
    write(rsp, file.path(out, 'rsp.txt'), 1, sep='\t')
    write(rownames(rsp), file.path(out, 'rsp.row'), 1)
    rdm <- c(rdm, "\nrsp.txt: N x 1 column vector of phenotypes, aligned with gmx\n")
    rdm <- c(rdm, "\trsp.row: row names of rsp, i.e., individual IDs.\n")

    ## true variance components
    vcs <- t(c(EPS=ret$eps, ret$vcs))
    write.table(vcs, file.path(out, 'vcs.txt'), sep='\t', quote=FALSE, row.names=FALSE)
    rdm <- c(rdm, "\nvcs.txt: row vector of true variance components.\n")

    ## fixed effect coefficients
    fix <- t(ret$fix)
    write.table(fix, file.path(out, 'fix.txt'), sep='\t', quote=FALSE, row.names=FALSE)
    rdm <- c(rdm, "\nfix.txt: row vector of true fixed effect coefficients, with intercept.\n")

    ## fixed effect matrix
    xmx <- ret$xmx
    write(t(xmx), file.path(out, 'xmx.txt'), ncol(xmx), sep='\t')
    write(rownames(xmx), file.path(out, 'xmx.row'), 1)
    write(colnames(xmx), file.path(out, 'xmx.col'), 1)
    rdm <- c(rdm, "\nxmx.txt: N x Q covariate matrix with intercept, aligned with gmx.\n")
    rdm <- c(rdm, "\txmx.row: row names of xmx, i.e., individual IDs.\n")
    rdm <- c(rdm, "\txmx.col: column names of xmx, i.e., the covariates.\n")

    ## variance component fitted by MINQUE
    if(!is.null(ret$fit))
    {
        fit <- t(ret$fit)
        write.table(fit, file.path(out, 'fit.txt'), sep='\t', quote=FALSE, row.names=FALSE)
        rdm <- c(rdm, "\nfit.mnq: row vector of MINQUE fitted variance components.\n")
    }

    ## default kernels
    if(!is.null(ret$knl))
    {
        . <- file.path(out, 'knl')
        dir.create(., FALSE, FALSE)
        knl <- with(ret, c(list(EPS=diag(nrow(gmx)), knl)))
        for(n in names(ret$knl))
        {
            k <- ret$knl[[n]]
            write(k, file.path(., paste0(n, '.txt')), ncol(k), sep='\t')
        }
        rdm <- c(rdm, "\nknl/: the working kernels to build models on,\n")
        rdm <- c(rdm, "\tEPS: diagnal/identity/noise kernel\n")
        rdm <- c(rdm, "\tLN1: standardized linear product kernel\n")
        rdm <- c(rdm, "\tLN2: standardized quadratic product kernel\n")
        rdm <- c(rdm, "\tLN3: standardized cubic product kernel\n")
        rdm <- c(rdm, "rows and columns of the kernels are aligned with gmx.\n")
    }

    ## write README
    writeLines(rdm, file.path(out, 'README'))
    
    invisible(NULL)
}

#' get simulated response
#'
#' @param dat genomic matrix, with training and testing masks.
#' @param frq frequency of functional variants
#' @param lnk link function to tranform the simulated signal;
#' @param eps size of noise
#' @param V the list of data generating kernels
#' @param het inter cohort heterogeneity
#' @param vc1 use these variance compoent for the shared effect,
#' the heterogeneity effect is still randomly generated.
#'
#' @return a list of lists, where the inner lists contains original genotypes and
#' generated response.
get.sim <- function(dat, ...)
{
    dot <- list(...)
    Q <- with(dat, length(unique(dvp[dvp > 0]))) # dvp groups
    R <- with(dat, length(unique(evl[evl > 0]))) # evl groups
    
    ## core effect
    dat <- within(dat,
    {
        if(sum(dvp > 0) > 0)
            dvp <- get.one(gmx[dvp > 0, ], ...)
        else
            rm(dvp)
        if(sum(evl > 0) > 0)
            evl <- get.one(gmx[evl > 0, ], ..., fmk=dvp$fmk)
        else
            rm(evl)
    })
    dat
}

sim <- function(fds, N=100, P=200, Q=1, R=1, ...)
{
    dot <- list(...)
    fit <- if(is.null(dot$fit)) FALSE else dot$fit

    ## get the list of R dataset
    rds <- get.rds(fds)
    
    ## the list of genomic data
    gno <- lapply(rds, readRDS)

    ## list of genomic matrices
    gmx <- get.gmx(gno, N, P, Q, R)
    
    ## generate simulation data
    sim <- get.sim(gmx, ...)

    ## fit a VCM by MINQUE?
    if(fit)
    {
        if(!is.null(sim$dvp))
            sim$dvp <- within(sim$dvp, fit <- mnq(rsp, knl, xmx, itr=40, tol=1e-5)$par)
        if(!is.null(sim$evl))
            sim$evl <- within(sim$evl, fit <- mnq(rsp, knl, xmx, itr=40, tol=1e-5)$par)
    }

    ## return
    if(!is.null(dot$out))
    {
        if(!is.null(sim$dvp))
            write.one(sim$dvp, paste0(dot$out, '.dvp'))
        if(!is.null(sim$evl))
            write.one(sim$evl, paste0(dot$out, '.evl'))
    }
    sim
}

main.gen <- function()
{
    ## ret <- lapply(1:10, function(i)
    ## {
    ##     fix <- round(c(X00=1, rnorm(2)), 1)
    ##     vcs <- round(rep(0.8, 3) + rchisq(3, 1) * 0.2, 1)
    ##     out <- sprintf("sim/gen/N03/%03d", i)
    ##     print(out)
    ##     sim('sim/dat', N=2000, P=4000, R=0, vcs=vcs, fix=fix, out=out, oks=~LN3, dks=~LN3, fit=TRUE)
    ## })

    ## ret <- lapply(1:10, function(i)
    ## {
    ##     fix <- round(c(X00=1, rnorm(2)), 1)
    ##     vcs <- round(rep(0.8, 2) + rchisq(2, 1) * 0.2, 1)
    ##     out <- sprintf("sim/gen/N02/%03d", i)
    ##     print(out)
    ##     sim('sim/dat', N=2000, P=4000, R=0, vcs=vcs, fix=fix, out=out, oks=~LN2, dks=~LN2, fit=TRUE)
    ## })

    ret <- lapply(1:10, function(i)
    {
        fix <- round(c(X00=1, rnorm(2)), 1)
        vcs <- round(1.0 * 0.8 + rchisq(1, 1) * 0.2, 1)
        out <- sprintf("sim/gen/N01/%03d", i)
        print(out)
        sim('sim/dat', N=2000, P=4000, R=0, vcs=vcs, fix=fix, out=out, oks=~LN1, dks=~LN1, fit=TRUE)
    })

    ret
}
