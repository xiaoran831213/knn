library(dplyr)
library(ggplot2)
source('R/hlp.R')

#' read training history from a directory
#' for simulation only
read.hist <- function(fdr, ..., par=FALSE, dvt=FALSE)
{
    hst <- list()
    rpt <- list()
    for(f in dir(c(fdr, ...), '.rds$', full.names=TRUE))
    {
        print(f)
        f <- readRDS(f)
        r <- f$rpt %>% as.tbl %>% select(N:bsz, mtd:usz)
        h <- f$hst %>% as.tbl %>% select(N:usz, bsz, i:par.phi, eve) %>%
            rename(luy=par.luy, lxu=par.lxu, phi=par.phi, val=eve)
        
        ## null error
        e <- r %>% filter(mtd=='nul', dat=='evl', key=='L2') %>% select(val) %>% unlist
        ## LMM error
        l <- r %>% filter(mtd=='lmm', dat=='evl', key=='L2') %>% select(val) %>% unlist

        ## relative error over null
        h <- h %>% mutate(rel=val/e, lmm=l/e)

        hst <- c(hst, list(h))
        rpt <- c(rpt, list(r))
    }
    hst <- bind_rows(hst)
    rpt <- bind_rows(rpt)
    
    ## separate configurations and recordings
    hst <- hst %>% group_by_at(vars(N:ep)) %>%
        summarise(sec=mean(sec), err=mean(rel), lmm=mean(lmm),
                  luy=mean(luy), lxu=mean(lxu), phi=mean(phi),
                  rep=n()) %>% ungroup
    hst <- hst %>% group_by_at(vars(N:usz)) %>% mutate(lmm=mean(lmm)) %>% ungroup
    rpt <- rpt %>% group_by_at(vars(-val)) %>% summarize(val=mean(val), rep=n()) %>% ungroup
    
    ## standardized epoch
    list(h=hst, r=rpt)
}

read.sim1 <- function(fdr, ..., par=FALSE, dvt=FALSE)
{
    hst <- list()
    rpt <- list()
    for(f in dir(c(fdr, ...), '.rds$', full.names=TRUE))
    {
        print(f)
        f <- readRDS(f)
        r <- f$rpt %>% as.tbl %>% select(-lr0:-max.itr)
        cyh.lmm <- r %>% filter(mtd=='lmm', dat=='dvp', key=='cyh') %>% select(val) %>% unlist
        mse.lmm <- r %>% filter(mtd=='lmm', dat=='dvp', key=='mse') %>% select(val) %>% unlist
        nlk.lmm <- r %>% filter(mtd=='lmm', dat=='dvp', key=='nlk') %>% select(val) %>% unlist
        cfg <- r %>% filter(row_number() == 1) %>% select(N:mtd)

        ## turn absolute performance into relativities to LMM
        PHI.kdn <- f$hst$par %>% lapply(`[[`, 'tuy') %>% sapply(`[`, 1, 1) %>% exp
        h <- cbind(cfg, f$hst$stt) %>% as.tbl %>%
            mutate(mnl=mnl/nlk.lmm, mse=mse/mse.lmm, cyh=cyh/cyh.lmm,
                   lxy=lxy/nlk.lmm, luy=luy/nlk.lmm, lxu=lxu/nlk.lmm, phi=PHI.kdn/PHI)

        hst <- cl(hst, h)
        rpt <- cl(rpt, r)
    }
    hst <- bind_rows(hst)
    rpt <- bind_rows(rpt)
    
    ## separate configurations and recordings
    hst.val <- hst %>% group_by_at(vars(N:i)) %>% summarise_all(mean) %>% ungroup
    hst.rep <- hst %>% group_by_at(vars(N:i)) %>% summarise(rep=n())
    hst <- inner_join(hst.val, hst.rep)

    rpt <- rpt %>% group_by_at(vars(-val)) %>% summarize(val=mean(val), rep=n()) %>% ungroup
    
    ## standardized epoch
    list(h=hst, r=rpt)
}

plot.sim1 <- function(aggr, out=NULL, xlim=NULL, ylim=NULL, smooth=FALSE, ...)
{
    rel <- aggr$h %>% filter(rep > max(rep) * .8)
    grp <- group_by_at(rel, vars(N:mtd))
    gnm <- length(group_size(grp))
    
    ## plotting device
    if(!is.null(out))
        png(out, width=1024 * gnm, height=1024)

    ## variables
    vn <- c('mse', 'cyh', 'phi', 'mnl', 'lxy', 'luy', 'lxu')
    gv <- group_vars(grp)

    ## shared title terms, and group specific title terms
    un <- sapply(rel[, gv], function(.) all(. == .[1]))
    tt <- paste0(gv[un], '=', unlist(rel[1, gv][un]), collapse=' ')
    gt <- gv[!un]
    
    lc <- rainbow(length(vn))           # line colors
    names(lc) <- vn
    lt <- c(1, 1, 1, 1, 2, 2, 2)        # line types
    names(lt) <- vn
    if(is.null(xlim))
        xlim <- with(rel, range(i))
    if(is.null(ylim))
        ylim <- c(0, 3)
    par(mfrow=c(1, gnm))
    .d <- do(grp,
    {
        plot(0, 0, 'n', xlim, ylim,
             xlab='training steps',
             ylab='performance relative to LMM')
        sbt <- paste0(gt, '=', unlist(.[1, gt]), collapse=' ')
        x <- .$i
        for(v in vn)
        {
            y <- unlist(.[v])
            if(lt[v] == 1 && smooth)
                y <- lowess(x, y, .05)[[2]]
            lines(x, y, col=lc[v], lty=lt[v], lwd=2)
        } 
        abline(1, 0)
        legend("topright", legend=vn, lwd=2, lty=lt, col=lc)
        title(sbt)
        .
    })
    if(!is.null(out))
        dev.off()
    NULL
}

#' plot training history
#' for simulation only
plot.hist.time <- function(agg, xlim=NULL, ylim=NULL)
{
    ## use simple strings to descript the kernels used
    agg <- agg %>% mutate(kcv=sprintf('%s(%s) ~ %s | %s', lnk, ycv, ikn, bkn)) %>%
        select(-lnk, -ycv, -bkn, -ikn)

    if(length(xlim) > 1)
        agg <- agg %>% filter(xlim[1] <= sec & sec <= xlim[2])
    else if(length(xlim) > 0)
        agg <- agg %>% filter(sec <= xlim)
    agg <- agg %>% mutate(PHI=as.factor(PHI), M=as.factor(M))
    
    g <- ggplot(agg)
    g <- g + geom_line(aes(x=sec, y=err, color=M))
    g <- g + geom_line(aes(x=sec, y=lmm, color=M))
    ## g <- g + geom_line(aes(x=sec, y=phi, color=M), linetype=2)
    ## g <- g + facet_grid(M~kcv)
    g <- g + facet_wrap(~PHI, nrow=1)
    if(!is.null(ylim))
        g <- g + ylim(ylim)
    ## else
    ##     g <- g + ylim(NA, 1.0)
    
    ## g <- g + scale_color_discrete(
    ## name='size of\nmini-batch:', breaks=c("128", "16"), labels = c("N/A", "16"))

    g <- g + xlab('Training Time (seconds)')
    g <- g + ylab('Relative Testing Error')
    g <- g + ggtitle('KDNN Performance by Training Time', 'the Effect of Noise:')
    g
}

plot.hist.pars <- function(agg, xlim=NULL, ylim=NULL)
{
    ## use simple strings to descript the kernels used
    agg <- agg %>% mutate(kcv=sprintf('%s(%s) ~ %s | %s', lnk, ycv, ikn, bkn)) %>%
        select(-lnk, -ycv, -bkn, -ikn)

    if(length(xlim) > 1)
        agg <- agg %>% filter(xlim[1] <= sec & sec <= xlim[2])
    else if(length(xlim) > 0)
        agg <- agg %>% filter(sec <= xlim)
    agg <- agg %>% mutate(nsi=as.factor(PHI), bsz=as.factor(bsz),
                          mmt=as.factor(mmt), usz=as.factor(usz), M=as.factor(M))
    
    g <- ggplot(agg)
    g <- g + geom_line(aes(x=sec, y=lxu+phi, color='LXU'))
    g <- g + geom_line(aes(x=sec, y=luy+phi, color='LUY'))
    g <- g + geom_line(aes(x=sec, y=phi-log(PHI), color='PHI'))
    g <- g + facet_grid(M~nsi)
    if(!is.null(ylim))
        g <- g + ylim(ylim)
    ## else
    ##     g <- g + ylim(NA, 1.0)
    
    ## g <- g + scale_color_discrete(
    ## name='size of\nmini-batch:', breaks=c("128", "16"), labels = c("N/A", "16"))

    g <- g + xlab('Training Time (seconds)')
    g <- g + ylab('exp(Parameters / Truth)')
    g <- g + ggtitle('KDN Consistency by Training Time', 'the Effect of Noise:')
    g
}
