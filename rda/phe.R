## phenotype handlers

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

## get design matrix
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
