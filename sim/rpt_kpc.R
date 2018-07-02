library(dplyr)
readRpt <- function(fdr, ..., ref=FALSE)
{
    rpt <- list()
    for(f in dir(c(fdr, ...), '.rds$', full.names=TRUE))
    {
        print(f)
        r <- try(readRDS(f))
        if(inherits(r, 'try-error'))
            next

        if(ref)
        {
            ## use null model as a measuring stick for MSE
            err <- subset(r, key %in% c('mse', 'nlk'))
            nul <- unlist(subset(err, mtd=='nul', val))
            err <- subset(err, mtd!='nul')
            err <- within(err, val <- val/nul)

            ## use GCTA as a measuring  stick for time elapsed.
            rtm <- subset(r, key=='rtm')
            gct <- unlist(subset(rtm, mtd=='gct', val))
            rtm <- subset(rtm, mtd!='gct')
            rtm <- within(rtm, val <- val/gct)

            r <- rbind(err, rtm)
        }
        rpt <- c(rpt, list(r))
    }
    rpt <- do.call(rbind, rpt)
    rpt
}

bxpRpt <- function(rpt, out=NULL, ...)
{
    ## use MSE to measure performance
    rpt <- subset(rpt, key=='mse', -c(key, seed))
    rpt <- subset(rpt, !mtd %in% c('nul', 'gmx', 'fmx'))

    ## simulation configurations
    grp <- subset(rpt, se=-c(mtd, val, dat))

    ## level 1 title shared by all subplot
    idx <- sapply(grp, function(x) length(unique(x)) == 1)
    hd1 <- unlist(grp[1, idx])
    hd1 <- paste0(names(hd1), '=', hd1, collapse=', ')

    ## split the data by configuration
    rpt <- split(rpt, grp)
    grp <- grp[, !idx, drop=FALSE]
    if(length(grp) < 1)
        grp <- list(grp)
    else
        grp <- lapply(split(grp, grp), unique)
    ngp <- length(rpt)
    
    ## plot
    if(!is.null(out))
        png(out, 4 * ngp, 4 * 2, units='in', res=300)

    par(mfcol=c(2, ngp), mar=c(4, 4, 1, 0) + .2, mgp=c(3, 1.3, 0), cex.lab=1.0, cex.axis=1.0)
    par(oma=c(0, 0, 2, 0))
    for(i in 1:ngp)
    {
        dat.dvp <- filter(rpt[[i]], dat=='dvp')
        dat.evl <- filter(rpt[[i]], dat=='evl')

        ## level 2 title specific to each subplot
        hd2 <- unlist(grp[[i]])
        hd2 <- paste0(names(hd2), '=', hd2, collapse=', ')
        xlb <- 'methods'

        boxplot(val~mtd, data=dat.dvp, ylab=ifelse(i==1, 'MSE of Fitting', ''), main=hd2, ...)
        boxplot(val~mtd, data=dat.evl, ylab=ifelse(i==1, 'MSE of Testing', ''), xlab=xlb, ...)
    }
    ## main title
    title(hd1, outer=TRUE)

    if(!is.null(out))
        dev.off()
    invisible(NULL)
}


#' write tab delimited file
WT <- function(d, f, ...) write.table(d, f, quote=FALSE, sep='\t', row.names=FALSE)

## boxplot of running time
bxpRtm <- function(rpt, out=NULL)
{
    ## get running time
    rpt <- subset(rpt, key=='rtm', -c(key, seed, dat))
    rpt <- subset(rpt, !mtd %in% c('nul', 'gmx', 'fmx'))
    ylm1 <- c(0, with(rpt, max(val)))
    
    ## simulation configurations
    grp <- subset(rpt, se=-c(mtd, val))

    ## level 1 title shared by all subplot
    idx <- sapply(grp, function(x) length(unique(x)) == 1)
    hd1 <- unlist(grp[1, idx])
    hd1 <- paste0(names(hd1), '=', hd1, collapse=', ')

    ## split the data by configuration
    rpt <- split(rpt, grp)
    grp <- grp[, !idx, drop=FALSE]
    if(length(grp) < 1)
        grp <- list(grp)
    else
        grp <- lapply(split(grp, grp), unique)
    ngp <- length(rpt)
    
    ## plot
    if(!is.null(out))
        png(out, 4 * ngp, 4 * 2, units='in', res=300)

    par(mfcol=c(2, ngp), mar=c(4, 4, 1, 0) + .2, mgp=c(3, 1.3, 0), cex.lab=1.0, cex.axis=1.0)
    par(oma=c(0, 0, 2, 0))
    for(i in 1:ngp)
    {
        dat <- rpt[[i]]

        ## level 2 title specific to each subplot
        hd2 <- unlist(grp[[i]])
        hd2 <- paste0(names(hd2), '=', hd2, collapse=', ')
        xlb <- 'methods'
        ylm2 <- c(0, 10)
        boxplot(val~mtd, data=dat, ylab=ifelse(i==1, 'Time Elapsed', ''), main=hd2, ylim=ylm1)
        abline(1, 0, col='red')
        boxplot(val~mtd, data=dat, ylab=ifelse(i==1, 'Time Elapsed', ''), xlab=xlb, ylim=ylm2)
        abline(1, 0, col='red')
    }
    ## main title
    title(hd1, outer=TRUE)

    if(!is.null(out))
        dev.off()
    invisible(NULL)
}

tm <- function()
{
    dir.create('rpt/2018_06_08')
    r1 <- readRpt('sim/s01'); 
    r2 <- readRpt('sim/s02'); 
    r3 <- readRpt('sim/s03'); 
    r4 <- readRpt('sim/s04'); 

    bxpRtm(r1, out='rpt/2018_06_08/t01_bxp.png')
    bxpRtm(r2, out='rpt/2018_06_08/t02_bxp.png')
    bxpRtm(r3, out='rpt/2018_06_08/t03_bxp.png')
    bxpRtm(r4, out='rpt/2018_06_08/t04_bxp.png')
}

er <- function()
{
    dir.create('rpt/2018_06_08')
    r1 <- readRpt('sim/s01');
    r2 <- readRpt('sim/s02');
    r3 <- readRpt('sim/s03');
    r4 <- readRpt('sim/s04');

    bxpRpt(r1, out='rpt/2018_06_08/e01_bxp.png', ylim=c(.2, 1.1))
    bxpRpt(r2, out='rpt/2018_06_08/e02_bxp.png', ylim=c(.0, 0.015))
    bxpRpt(r3, out='rpt/2018_06_08/e03_bxp.png', ylim=c(.25, 0.55))
    bxpRpt(r4, out='rpt/2018_06_08/e04_bxp.png', ylim=c(.25, 0.55))
}
