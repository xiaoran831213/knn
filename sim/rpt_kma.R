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

    ## mq1 <- subset(rpt, mtd == 'mq1')
    ## mq2 <- subset(rpt, mtd == 'mq2')
    ## mq3 <- subset(rpt, mtd == 'mq3')
    ## mq4 <- subset(rpt, mtd == 'mq4')
    ## gct <- subset(rpt, mtd == 'gct')

    ## mq1 <- within(mq1, {mtd <- 'MNQ-M1'})
    ## mq2 <- within(mq2, {mtd <- 'MNQ-M2'})
    ## mq3 <- within(mq3, {mtd <- 'MNQ-BT'})
    ## mq4 <- within(mq4, {mtd <- 'MNQ-WH'})

    ## mnq <- rbind(mq1, mq2, mq3, mq4)
    ## gct <- within(gct, {mtd <- 'GCTA'})
    ## rpt <- rbind(mnq, gct)

    ## ylim
    ylm <- with(rpt, c(min(val), min(2, max(val))))
    
    ## simulation configurations
    grp <- subset(rpt, se=c(ejt))
    
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
        png(out, ngp * 4, 1 * 4, units='in', res=300)

    par(mfcol=c(1, ngp), mar=c(2.25, 2.35, 1.2, 0), mgp=c(1.5, 0.5, 0), cex.lab=1.0, cex.axis=1.0)
    par(oma=c(0, 0, 2, 0) +.1)
    for(i in 1:ngp)
    {
        dat <- rpt[[i]]

        ## level 2 title specific to each subplot
        hd2 <- unlist(grp[[i]])
        hd2 <- paste0(names(hd2), '=', hd2, collapse=', ')
        xlb <- 'methods'
        ylb <- 'MSE'

        boxplot(val~mtd, dat, xlab=xlb, ylab=ylb, ylim=ylm, ...)
        abline(1, 0, col='red')
        title(hd2)
    }
    ## main title
    title(hd1, outer=TRUE, adj=0)

    if(!is.null(out))
        dev.off()
    invisible(NULL)
}

sm <- function()
{
    dir.create('rpt/2018_06_20')
    s0 <- readRpt('sim/s00', ref=1);
    s1 <- readRpt('sim/s01', ref=1);
    s2 <- readRpt('sim/s02', ref=1);
    s3 <- readRpt('sim/s03', ref=1);

    bxpRpt(s0, out='rpt/2018_06_20/s00_mse.png')
    bxpRpt(s1, out='rpt/2018_06_20/s01_mse.png')
    bxpRpt(s2, out='rpt/2018_06_20/s02_mse.png')
    bxpRpt(s3, out='rpt/2018_06_20/s03_mse.png')

    bxpRtm(s0, out='rpt/2018_06_20/s00_rtm.png')
    bxpRtm(s1, out='rpt/2018_06_20/s01_rtm.png')
    bxpRtm(s2, out='rpt/2018_06_20/s02_rtm.png')
    bxpRtm(s3, out='rpt/2018_06_20/s03_rtm.png')
}
