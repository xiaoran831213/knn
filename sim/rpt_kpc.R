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
    rpt <- rpt %>% filter(key=='mse') %>% select(-seed, -key)
    agg <- rpt %>% group_by_at(vars(N:key)) %>% summarize(val=mean(val)) %>% ungroup
    agg
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

## only on testing data
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


p1 <- function()
{
    rpt <- readRpt('sim/s01', ref.nul=TRUE)
    rpt <- subset(rpt, mtd %in% c("gct", "rop", "mq1", "mq2"))

    rpt <- within(rpt,
    {
        lnk <- "I"
        mtd[mtd=="gct"] <- "GCTA\nREML"
        mtd[mtd=="rop"] <- "R\nML"
        mtd[mtd=="mq1"] <- "MINQUE\norder1"
        mtd[mtd=="mq2"] <- "MINQUE\norder2"
        yks <- "Id, G'G"
        oks <- "Id, Gaussian(X)"
    })

    rpt <- subset(rpt, frq==0.1)
    bxpRpt(rpt, out='rpt/2018_05_31/s01_bxp.png')
    invisible(rpt)
}

p2 <- function()
{
    rpt <- readRpt('sim/s02', ref.nul=TRUE)
    rpt <- subset(rpt, mtd %in% c("gct", "rop", "mq1", "mq2"))

    rpt <- within(rpt,
    {
        lnk <- "sigmoid"
        mtd[mtd=="gct"] <- "GCTA\nREML"
        mtd[mtd=="rop"] <- "R\nML"
        mtd[mtd=="mq1"] <- "MINQUE\norder1"
        mtd[mtd=="mq2"] <- "MINQUE\norder2"
        yks <- "Id, G'G"
        oks <- "Id, X'X"
    })
    bxpRpt(rpt, out='rpt/2018_05_31/s02_bxp.png')
    invisible(rpt)
}


p3 <- function()
{
    rpt <- readRpt('sim/s03', ref.nul=TRUE)
    rpt <- subset(rpt, mtd %in% c("gct", "rop", "mq1", "mq2"))

    rpt <- within(rpt,
    {
        lnk <- "sin(1 * pi * x)"
        mtd[mtd=="gct"] <- "GCTA\nREML"
        mtd[mtd=="rop"] <- "R\nML"
        mtd[mtd=="mq1"] <- "MINQUE\norder1"
        mtd[mtd=="mq2"] <- "MINQUE\norder2"
        yks <- "Id, G'G"
        oks <- "Id, X'X"
    })
    bxpRpt(rpt, out='rpt/2018_05_31/s03_bxp.png')
    invisible(rpt)
}

p4 <- function()
{
    rpt <- readRpt('sim/s04', ref.nul=TRUE)
    rpt <- subset(rpt, mtd %in% c("gct", "rop", "mq1", "mq2"))

    rpt <- within(rpt,
    {
        lnk <- "sin(2 * pi * x)"
        mtd[mtd=="gct"] <- "GCTA\nREML"
        mtd[mtd=="rop"] <- "R\nML"
        mtd[mtd=="mq1"] <- "MINQUE\norder1"
        mtd[mtd=="mq2"] <- "MINQUE\norder2"
        yks <- "Id, G'G"
        oks <- "Id, X'X"
    })
    bxpRpt(rpt, out='rpt/2018_05_31/s04_bxp.png')
    invisible(rpt)
}

#' write tab delimited file
WT <- function(d, f, ...) write.table(d, f, quote=FALSE, sep='\t', row.names=FALSE)


a1 <- function()
{
    r <- readRpt('sim/s01', ref.nul=TRUE)
    r <- r %>% filter(eps==.2, key=='mse', mtd!='fmx', mtd!='gmx') %>% select(-seed, -key)
    a <- r %>% group_by_at(vars(N:dat)) %>% summarize(val=mean(val)) %>% ungroup
    WT(a, 'rpt/2018_05_31/s01_agg.txt')
    a
}

a2 <- function()
{
    r <- readRpt('sim/s02', ref.nul=TRUE)
    r <- r %>% filter(eps==.2, key=='mse', mtd!='fmx', mtd!='gmx') %>% select(-seed, -key)
    a <- r %>% group_by_at(vars(N:dat)) %>% summarize(val=mean(val)) %>% ungroup
    WT(a, 'rpt/2018_05_31/s02_agg.txt')
    a
}

a3 <- function()
{
    r <- readRpt('sim/s03', ref.nul=TRUE)
    r <- r %>% filter(eps==.2, key=='mse', mtd!='fmx', mtd!='gmx') %>% select(-seed, -key)
    a <- r %>% group_by_at(vars(N:dat)) %>% summarize(val=mean(val)) %>% ungroup
    WT(a, 'rpt/2018_05_31/s03_agg.txt')
    a
}

a4 <- function()
{
    r <- readRpt('sim/s04', ref.nul=TRUE)
    r <- r %>% filter(eps==.2, key=='mse', mtd!='fmx', mtd!='gmx') %>% select(-seed, -key)
    a <- r %>% group_by_at(vars(N:dat)) %>% summarize(val=mean(val)) %>% ungroup
    WT(a, 'rpt/2018_05_31/s04_agg.txt')
    a
}

## boxplot of running time
bxpRtm <- function(rpt, out=NULL)
{
    ## get running time
    rpt <- subset(rpt, key=='rtm', -c(key, seed, dat))
    rpt <- subset(rpt, !mtd %in% c('nul', 'gmx', 'fmx'))
    
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
        ylm1 <- c(0, 50)
        ylm2 <- c(0, 10)
        boxplot(val~mtd, data=dat, ylab=ifelse(i==1, 'Time Elapsed', ''), main=hd2, ylim=ylm1)
        boxplot(val~mtd, data=dat, ylab=ifelse(i==1, 'Time Elapsed', ''), xlab=xlb, ylim=ylm2)
    }
    ## main title
    title(hd1, outer=TRUE)

    if(!is.null(out))
        dev.off()
    invisible(NULL)
}

tm <- function()
{
    r1 <- readRpt('sim/s01'); bxpRtm(r1, out='rpt/2018_05_31/t01_bxp.png')
    r2 <- readRpt('sim/s02'); bxpRtm(r2, out='rpt/2018_05_31/t02_bxp.png')
    r3 <- readRpt('sim/s03'); bxpRtm(r3, out='rpt/2018_05_31/t03_bxp.png')
    r4 <- readRpt('sim/s04'); bxpRtm(r4, out='rpt/2018_05_31/t04_bxp.png')
}
