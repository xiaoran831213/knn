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

    mq1 <- subset(rpt, mtd == 'mq1')
    mq2 <- subset(rpt, mtd == 'mq2')
    mq3 <- subset(rpt, mtd == 'mq3')
    mq4 <- subset(rpt, mtd == 'mq4')
    gct <- subset(rpt, mtd == 'gct')

    mq1 <- within(mq1, {mtd <- 'L-MNQ'; bsz <- N})
    mq2 <- within(mq2, {mtd <- 'L-MNQ'})
    mq3 <- within(mq3, {mtd <- 'P-MNQ'; bsz <- N})
    mq4 <- within(mq4, {mtd <- 'P-MNQ'})
    
    mnq <- rbind(mq1, mq2, mq3, mq4)
    gct <- within(gct, {mtd <- 'GCTA'; bsz <- N})
    rpt <- rbind(mnq, gct)
    
    ## ylim
    ylm <- with(rpt, c(min(val), min(2, max(val))))
    
    ## simulation configurations
    grp <- subset(rpt, se=-c(bsz, val, dat))

    ## level 1 title shared by all subplot
    idx <- sapply(grp, function(x) length(unique(x)) == 1)
    hd1 <- unlist(grp[1, idx])
    hd1 <- paste0(names(hd1), '=', hd1, collapse=', ')
    hd1 <- paste0('Mean square error relative to null model\n', hd1)

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

    par(mfcol=c(2, ngp), mar=c(2.25, 2.35, 1.2, 0), mgp=c(1.5, 0.5, 0), cex.lab=1.0, cex.axis=1.0)
    par(oma=c(0, 0, 2, 0) +.1)
    for(i in 1:ngp)
    {
        dat.dvp <- filter(rpt[[i]], dat=='dvp')
        dat.evl <- filter(rpt[[i]], dat=='evl')

        ## level 2 title specific to each subplot
        hd2 <- unlist(grp[[i]])
        hd2 <- paste0(names(hd2), '=', hd2, collapse=', ')
        xlb <- 'batch size'

        boxplot(val~bsz, data=dat.dvp,
                ylab=ifelse(i==1, 'MSE of Fitting', ''), main=hd2, ylim=ylm, ...)
        abline(1, 0, col='red')
        boxplot(val~bsz, data=dat.evl,
                ylab=ifelse(i==1, 'MSE of Testing', ''), xlab=xlb, ylim=ylm, ...)
        abline(1, 0, col='red')
    }
    ## main title
    title(hd1, outer=TRUE, adj=0)

    if(!is.null(out))
        dev.off()
    invisible(NULL)
}


#' write tab delimited file
WT <- function(d, f, ...) write.table(d, f, quote=FALSE, sep='\t', row.names=FALSE)

## boxplot of running time
bxpRtm <- function(rpt, out=NULL)
{
    ## use MSE to measure performance
    rpt <- subset(rpt, key=='rtm' & dat=='dvp', -c(key, seed, dat))

    mq1 <- subset(rpt, mtd == 'mq1')
    mq2 <- subset(rpt, mtd == 'mq2')
    mq3 <- subset(rpt, mtd == 'mq3')
    mq4 <- subset(rpt, mtd == 'mq4')

    mq1 <- within(mq1, {mtd <- 'L-MNQ'; bsz <- N})
    mq2 <- within(mq2, {mtd <- 'L-MNQ'})
    mq3 <- within(mq3, {mtd <- 'P-MNQ'; bsz <- N})
    mq4 <- within(mq4, {mtd <- 'P-MNQ'})
    
    rpt <- rbind(mq1, mq2, mq3, mq4)
    ylm1 <- unname(quantile(rpt$val, c(.01, .99)))
    
    ## simulation configurations
    grp <- subset(rpt, se=-c(bsz, val))

    ## level 1 title shared by all subplot
    idx <- sapply(grp, function(x) length(unique(x)) == 1)
    hd1 <- unlist(grp[1, idx])
    hd1 <- paste0(names(hd1), '=', hd1, collapse=', ')
    hd1 <- paste0('Training Time relative to GCTA\n', hd1)

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

    par(mfcol=c(2, ngp), mar=c(2.25, 2.25, 1.2, 0), mgp=c(1.5, 0.5, 0), cex.lab=1.0, cex.axis=1.0)
    par(oma=c(0, 0, 2, 0) +.1)
    for(i in 1:ngp)
    {
        dat <- rpt[[i]]

        ## level 2 title specific to each subplot
        hd2 <- unlist(grp[[i]])
        hd2 <- paste0(names(hd2), '=', hd2, collapse=', ')
        xlb <- 'batch size'
        ylm2 <- c(0, 1.1)
        boxplot(val~bsz, data=dat, ylab=ifelse(i==1, 'full range', ''), main=hd2, ylim=ylm1)
        abline(1, 0, col='red')
        boxplot(val~bsz, data=dat, ylab=ifelse(i==1, 'in [0, 1]', ''), xlab=xlb, ylim=ylm2)
        abline(1, 0, col='red')
    }
    ## main title
    title(main=hd1, outer=TRUE, adj=0)

    if(!is.null(out))
        dev.off()
    invisible(NULL)
}


bt <- function()
{
    dir.create('rpt/2018_06_08')
    b1 <- readRpt('sim/b01', ref=1);
    b2 <- readRpt('sim/b02', ref=1);
    b3 <- readRpt('sim/b03', ref=1);
    b4 <- readRpt('sim/b04', ref=1);

    bxpRpt(b1, out='rpt/2018_06_08/b01_mse.png')
    bxpRpt(b2, out='rpt/2018_06_08/b02_mse.png')
    bxpRpt(b3, out='rpt/2018_06_08/b03_mse.png')
    bxpRpt(b4, out='rpt/2018_06_08/b04_mse.png')

    bxpRtm(b1, out='rpt/2018_06_08/b01_rtm.png')
    bxpRtm(b2, out='rpt/2018_06_08/b02_rtm.png')
    bxpRtm(b3, out='rpt/2018_06_08/b03_rtm.png')
    bxpRtm(b4, out='rpt/2018_06_08/b04_rtm.png')
}

sm <- function()
{
    dir.create('rpt/2018_06_14')
    s0 <- readRpt('sim/s00', ref=1);
    s1 <- readRpt('sim/s01', ref=1);
    s2 <- readRpt('sim/s02', ref=1);
    s3 <- readRpt('sim/s03', ref=1);

    bxpRpt(s0, out='rpt/2018_06_14/s00_mse.png')
    bxpRpt(s1, out='rpt/2018_06_14/s01_mse.png')
    bxpRpt(s2, out='rpt/2018_06_14/s02_mse.png')
    bxpRpt(s3, out='rpt/2018_06_14/s03_mse.png')

    bxpRtm(s0, out='rpt/2018_06_14/s00_rtm.png')
    bxpRtm(s1, out='rpt/2018_06_14/s01_rtm.png')
    bxpRtm(s2, out='rpt/2018_06_14/s02_rtm.png')
    bxpRtm(s3, out='rpt/2018_06_14/s03_rtm.png')
}
