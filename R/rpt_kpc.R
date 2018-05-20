library(dplyr)
readRpt <- function(fdr, ..., ref.nul=FALSE)
{
    rpt <- list()
    for(f in dir(c(fdr, ...), '.rds$', full.names=TRUE))
    {
        print(f)
        r <- try(readRDS(f))
        if(inherits(r, 'try-error'))
            next

        if(ref.nul)
        {
            r <- subset(r, key!='rtm' & key!='cyh')
            n <- unlist(subset(r, mtd=='nul', val))
            r <- subset(r, mtd!='nul')
            r <- within(r, val <- val/n)
        }
        rpt <- c(rpt, list(r))
    }
    rpt <- do.call(rbind, rpt)
    rpt
}

aggRpt <- function(rpt)
{
    agg <- rpt %>% group_by_at(vars(N:key)) %>% summarize(val=mean(val)) %>% ungroup
}

bxpRpt <- function(rpt, out=NULL)
{
    ## use MSE to measure performance
    rpt <- subset(rpt, key=='mse', -c(key, seed))
    
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

    par(mfcol=c(2, ngp), mar=c(4, 4, 1, 0) + .2, cex.lab=1.2, cex.axis=1.2)
    par(oma=c(0, 0, 2, 0))
    for(i in 1:ngp)
    {
        dat.dvp <- filter(rpt[[i]], dat=='dvp')
        dat.evl <- filter(rpt[[i]], dat=='evl')

        ## level 2 title specific to each subplot
        hd2 <- unlist(grp[[i]])
        hd2 <- paste0(names(hd2), '=', hd2, collapse=', ')
        xlb <- 'methods'

        ylm <- c(0, 1.2)
        boxplot(val~mtd, data=dat.dvp, ylim=ylm, ylab=ifelse(i==1, 'MSE of Fitting', ''), main=hd2)
        boxplot(val~mtd, data=dat.evl, ylim=ylm, ylab=ifelse(i==1, 'MSE of Testing', ''), xlab=xlb)
    }
    ## main title
    title(hd1, outer=TRUE)

    if(!is.null(out))
        dev.off()
    invisible(NULL)
}

bxpRp1 <- function(rpt, out=NULL)
{
    ## use MSE to measure performance
    rpt <- subset(rpt, key=='mse' & dat=='evl', -c(key, dat, seed))
    
    ## plot
    if(!is.null(out))
        png(out, 4, 4, units='in', res=300)

    par(mar=c(4, 4, 1, 0) + .2, cex.lab=1, cex.axis=1, cex.main=.9)
    par(oma=c(0, 0, 2, 0))

    ## level 2 title specific to each subplot
    hd1 <- rpt[1, c('N', 'P', 'frq', 'lnk', 'eps')]
    hd2 <- rpt[1, c('yks', 'oks')]
    hd1 <- paste0(names(hd1), '=', hd1, collapse=', ')
    hd2 <- paste0(names(hd2), '=', hd2, collapse=', ')
    hdr <- paste0(hd1, '\n', hd2)
    xlb <- 'methods'
    
    ylm <- c(0, 1.2)
    boxplot(val~mtd, data=rpt, xlab=xlb, ylim=ylm, ylab='MSE of Testing')

    ## main title
    title(hdr, outer=TRUE)

    if(!is.null(out))
        dev.off()
    invisible(NULL)
}


main <- function()
{
    rpt <- readRpt('sim/s02', ref.nul=TRUE)
    rpt <- subset(rpt, mtd %in% c("gct", "rop", "mnq"))
    rpt <- within(rpt,
    {
        lnk <- "I"
        mtd[mtd=="gct"] <- "GCTA-REML"
        mtd[mtd=="rop"] <- "R-ML"
        mtd[mtd=="mnq"] <- "MINQUE"
        yks <- "e*I + s^2*G'G"
        oks <- "e*I + s^2*X'X"
    })

    rpt <- subset(rpt, frq==0.1)
    bxpRp1(rpt)
    ## bxpRpt(rpt, out='data/rpt/sm3_bxp.png')
    rpt
}
