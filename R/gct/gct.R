## GCTA wrappers

## GCTA REML
## y: response
## g: GRMs
## y and GRMS are assumed to be aligned.
## make sure the name 'gcta64' points to the executable
gcta.reml <- function(y, K, qcvr=NULL, dcvr=NULL, maxit=NULL)
{
    ## temporary directory
    tpd <- paste0("gcta", tempfile("", ''))
    if(!dir.create(tpd, FALSE, TRUE))
        stop('failed to create directory', GRM.dir)

    ## save GRM
    GRM <- K
    if(is.null(names(GRM)))
        names(GRM) <- sprintf('G%02d', seq(length(GRM)))
    K <- c(list(phi=diag(NROW(y))), K)

    GRM.dir <- file.path(tpd, 'grm')
    if(!dir.create(GRM.dir))
        stop('failed to create directory', GRM.dir)
    GRM.path <- file.path(GRM.dir, names(GRM))
    for(i in seq_along(GRM))
    {
        saveGRM(GRM.path[i], GRM[[i]])
        PF("save GRM: %d %16s %8d %8d\n", i, GRM.path[i], nrow(GRM[i]), ncol(GRM[i]))
    }
    mgrm.path <- file.path(tpd, 'grm.lst')
    write(GRM.path, mgrm.path)

    ## write phenotype and covariate
    id <- makeID(GRM[[1]])
    phe <- cbind(id, y=y)
    phe.path <- file.path(tpd, 'phe.txt')
    write.table(phe, phe.path, quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)

    ## compose GCTA command
    exe <- 'gcta64'                     # binary
    fun <- '--reml'                     # REML
    if(length(GRM) == 1)                # GRM(s)
        mgm <- paste('--grm', GRM.path[[1]])
    else
        mgm <- paste('--mgrm', mgrm.path)
    phe <- paste('--pheno', phe.path)
    out <- paste('--out', file.path(tpd, 'out'))

    opt <- paste('--reml-pred-rand', '--reml-no-lrt')
    if(!is.null(maxit))
        opt <- paste(opt, '--reml-maxit', maxit)
    cmd <- paste(exe, fun, mgm, phe, out, opt)

    ## execute command
    ext <- system(cmd, show.output.on.console=FALSE)
    if(ext != 0)
        stop(cmd, ' GCTA existed with non-zero.')

    ## parse the output
    out <- gcta.parse(tpd)
    
    ## summary
    h <- rowSums(out$blp[-1:-3])
    v <- cmb(K, out$vcs$var)[[1]]
    ## mse <- mean((y - h)^2)              # mean squre error
    mse <- mean(out$blp[, 3]^2)         # estimate residual
    cyh <- cor(y, h)
    nlk <- nlk(y, v)
    rpt <- DF(key=c('mse', 'nlk', 'cyh', 'rtm'), val=c(mse, nlk, cyh, out$rtm))

    ## remove temporary directory
    unlink(tpd, TRUE, TRUE)

    ## pack up and return
    c(list(cmd=cmd, ext=ext, rpt=rpt), out)
}

gcta.parse <- function(tpd)
{
    ## --- parse *.hsq for variance ---
    hsq.path <- file.path(tpd, 'out.hsq')
    hsq <- readLines(hsq.path)

    ## variance components, and their names
    reg <- "(^Source|^V[(])"
    vcs <- grep(reg, hsq, value=TRUE)
    hsq <- grep(reg, hsq, value=TRUE, invert=TRUE)
    vcs <- read.delim(text=vcs)
    names(vcs) <- c('par', 'var', 'se')
    vcs <- within(vcs, par <- sub("V[(](.*)[)]", "\\1", par))
    vcs <- subset(vcs, !grepl('Vp$', par))
    vcs <- vcs[c(nrow(vcs), 1:(nrow(vcs) - 1)), ]
    vcs$par[1] <- 'PHI'
    
    ## --- parse *.indi.blp for predictions on training data ---
    blp.path <- file.path(tpd, 'out.indi.blp')
    blp <- read.table(blp.path)
    blp <- blp[, c(1, 2, seq(4, ncol(blp), 2))] # effects
    blp <- blp[, c(1:2, ncol(blp), seq(3, ncol(blp) - 1))]
    names(blp)[-(1:2)] <- vcs$par

    ## --- parse *.log for time elapsed
    log.path <- file.path(tpd, 'out.log')
    log <- readLines(log.path)
    rtm <- grep("^Computational time", log, value=TRUE)
    rtm <- as.numeric(sub("^.*: ([0-9.]*) .*$", "\\1", rtm))

    ## return
    par <- vcs$var
    list(vcs=vcs, hsq=hsq, blp=blp, par=par, log=log, rtm=rtm)
}
