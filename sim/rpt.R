library(ggplot2)
library(dplyr)
source('R/hlp.R')

plt.cor <- function(dat)
{
    d1 <- list()
    nm <- names(dat)
    for(i in seq_along(dat))
    {
        for(j in seq_along(dat))
        {
            p1 <- data.frame(dx=dat[, i], dy=dat[, j], vx=nm[i], vy=nm[j])
            d1 <- c(d1, list(p1))
        }
    }
    d1 <- do.call(rbind, d1)

    gc <- ggplot(d1, aes(x=dx, y=dy))
    gc <- gc + geom_point(size=.5)
    gc <- gc + facet_grid(vx~vy)

    th <- theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    gc <- gc + th
    gc
}

readRPT <- function(..., ref=1)
{
    rpt <- list()
    dot <- c(...)
    for(f in dir(dot, '.rds$', full=TRUE))
    {
        print(f)
        r <- try(readRDS(f))
        if(inherits(r, 'try-error'))
            next

        evl <- subset(r, dat == 'evl')
        dvp <- subset(r, dat == 'dvp')
        bia <- subset(r, dat == 'bia')
        
        ## MSE
        mse <- subset(evl, key == 'mse')
        if(ref)
        {
            nul <- subset(mse, mtd == 'NUL')$val
            mse <- within(mse, val <- val / nul)
        }

        ## likelihood
        nlk <- subset(evl, key == 'nlk')
        if(ref)
        {
            nul <- subset(nlk, mtd == 'NUL')$val
            ## ldt <- within(subset(evl, key == 'ldt'), val <- exp(val - nul))
            ## yay <- within(subset(evl, key == 'yay'), val <- exp(val - nul))
            nlk <- within(nlk, val <- exp(2 * (val - nul)))
        }

        ## bias of epsilon
        
        bs0 <- within(subset(bia, key == 'eps'), {key <- 'bs0'; dat='evl'})
        if(ref)
            bs0 <- within(bs0, val  <- val / 4 + 0.5)

        ## R^2
        rsq <- subset(evl, key == 'rsq')

        ## running time
        rtm <- subset(dvp, key == 'rtm')
        if(ref)
        {
            rtm <- within(rtm, val <- val / median(val) / 2)
        }
        rpt <- c(rpt, list(rbind(bia, rtm, mse, nlk, bs0, rsq)))
    }
    rpt <- do.call(rbind, rpt)
    
    rpt <- within(rpt, {H <- N * R; N <- N * Q})
    rpt <- within(rpt, rm(R, Q, seed))
    rpt <- within(rpt,
    {
        tag <- as.factor(tag)
        if('ref' %in% levels(tag))
           tag <- relevel(tag, 'ref')
    })
    rpt
}

## infer title from repeatative columns.
ttl <- function(d)
{
    i <- sapply(d, function(x) length(unique(x)) == 1)
    u <- as.list(d[1, i, drop=FALSE])
    paste0(names(u), '=', u, collapse=', ')
}

## fetch reports
d0 <- function(name, use.cache=TRUE, ...)
{
    fdr <- file.path('sim/run', name)
    rds <- paste0(fdr, '.rds')
    if(file.exists(rds) && use.cache)
        dat <- readRDS(rds)
    else
    {
        dat <- readRPT(fdr, ...)
        saveRDS(dat, rds)
    }
    print(rds)
    dat
}

## plot reports
prpt <- function(name, cache=TRUE, errs=TRUE, bias=TRUE, ...)
{
    d <- as_tibble(d0(name, cache, ...))
    d <- filter(d, !mtd %in% c('NUL'))
    r <- list()

    th <- theme(
        axis.title.x=element_blank(), axis.title.y=element_blank(), 
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="red", fill="#CCCCFF"))
    ## report: errors
    if(errs)
    {
        rpt <- filter(d, key %in% c('bs0', 'rsq', 'mse', 'nlk', 'rtm'))
        rpt <- mutate(rpt, val=pmin(val, 1.00), val=pmax(val, 0.00))
        g <- ggplot(rpt, aes(x=mtd, y=val))
        g <- g + geom_boxplot()
        g <- g + facet_grid(tag~key, scales='free')
        ## g <- g + ylim(c(0, 1))
        g <- g + ggtitle(ttl(rpt))
        g <- g + th
        f <- paste0(file.path('~/img', name), '_rpt.png'); print(f)
        h <- length(unique(rpt$tag))
        w <- length(unique(rpt$key))
        ggsave(f, g, width=1.9 * w + .2, height=1 * h + .2, scale=1.3)
        r <- c(r, rpt=g)
    }

    if(bias)
    {
        ## report: bias
        bia <- filter(d, dat=='bia')
        bia <- mutate(bia, key=sub('^.*[.]', '', key))
        bia <- mutate(bia, val=pmin(val, 1.99), val=pmax(val, -1.99))
        g <- ggplot(bia, aes(x=key, y=val))
        g <- g + geom_boxplot()
        g <- g + facet_grid(tag~mtd)
        ## g <- g + ylim(c(-2, 2))
        g <- g + ggtitle(ttl(bia))
        g <- g + th
        w <- length(unique(rpt$mtd))
        h <- length(unique(rpt$tag))
        f <- paste0(file.path('~/img', name), '_bia.png'); print(f)
        ggsave(f, g, width=1.9 * w + .2, height=1 * h + .2, scale=1.3)
        r <- c(r, bia=g)
    }

    ## return
    invisible(r)
}

## labeller
lb1 <- function (labels, multi_line = TRUE) 
{
    labels <- label_value(labels, multi_line = multi_line)
    dc <- c(
        `ref`="bold(y)*'='*bold(alpha)+bold(epsilon)",
        `bin`="binomial",
        `chi`="'chi-square'",
        `poi`="Poisson",
        `exp`="exponential",
        `hy1`="bold(y)*'='*over('|'*bold(alpha)*'|', 1 + '|'*bold(alpha)*'|') + bold(epsilon)",
        `rc1`="bold(y)*'='*bold('|'*alpha*'|')*e^bold('|'*alpha*'|') + bold(epsilon)",
        `e^1`="bold(y)*'='*bold(alpha) + bold(epsilon)",
        `e^2`="bold(y)*'='*bold(alpha)^2 + bold(epsilon)",
        `e^3`="bold(y)*'='*bold(alpha)^3 + bold(epsilon)",
        `g:2`="'2-way interaction'",
        `g:3`="'3-way interaction'",
        `g*2`="'2-way + quadratic'",
        `g*3`="'3-way + cubic'",
        `bs0`="hat(symbol(sigma))[0]^2 - symbol(sigma)[0]^2",
        `mse`="MSE(bold(y), hat(bold(y)))",
        `rsq`="COR(bold(y), hat(bold(y)))^2",
        `rtm`="time ~~ (sec)")
    if (multi_line)
    {
        lapply(unname(labels), lapply, function(values)
        {
            if(values %in% names(dc))
                values <- dc[[values]]
            c(parse(text = as.character(values)))
        })
    }
    else
    {
        lapply(labels, function(values)
        {
            values <- paste0("list(", values, ")")
            lapply(values, function(expr) c(parse(text = expr)))
        })
    }
}

pabs <- function(d, o=NULL, fill=TRUE, bat=FALSE)
{
    th <- theme(
        axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), 
        strip.text.x = element_text(size=12, face="bold"),
        strip.text.y = element_text(size=12, face="bold"),
        strip.background = element_rect(colour="red", fill="#CCCCFF"))
    if(fill)
    {
        th <- th + theme(legend.title=element_blank(), legend.position='bottom')
        th <- th + theme(axis.text.x=element_blank())
    }
    
    ## method dictionary
    if(!bat)
        d <- subset(d, !grepl('^B', mtd))
    d <- within(d, mds <- c(
                  GCT="GCTA",
                  MNQ="KNN",
                  BMQ="KNN-BAT")[mtd])
    d <- within(d, mds <- as.factor(mds))
    d <- within(d, relevel(mds, 'GCTA'))
    
    ## report: errors grouped by (tag, key)
    d <- filter(d, key %in% c('bs0', 'rsq', 'mse', 'rtm'))

    ## cap the outliers
    d <- d %>% group_by(tag, key) %>%
        mutate(val=pmin(val, quantile(val, .975))) %>%
        mutate(val=pmax(val, quantile(val, .025))) %>% ungroup

    ## dimensional of sub-plots
    ntag <- length(unique(d$tag))
    nkey <- length(unique(d$key))

    ## plot
    if(fill)
        g <- ggplot(d, aes(x=mtd, y=val, fill=mds))
    else
        g <- ggplot(d, aes(x=mtd, y=val))
    g <- g + geom_boxplot()
    ## if(bat)
    ##     g <- g + geom_boxplot(data=subset(d, grepl('^BM.', mtd)))
    ## else
    ##     g <- g + geom_boxplot(data=subset(d, grepl('^BM.', mtd)), alpha=0, linetype=0)
    g <- g + facet_grid(key ~ tag, scales='free_y', labeller=lb1)
    g <- g + th
    g <- g + ggtitle(ttl(d))
    ## g <- g + scale_fill_manual(breaks=c('GCTA', 'KNN', 'KNN-BAT'),
    ##                            values=c("red", "blue", "green"))

    ## output
    if(!is.null(o))
    {
        o <- paste0(file.path('~/img', o), '.png')
        print(o)
        ggsave(o, g, width=1.9 * ntag + .2, height=1.0 * nkey + .2, scale=1.5)
    }
    ## return
    invisible(d)
}


rpt1 <- function(cache=FALSE)
{
    . <- subset(d0('bi0', cache), mtd != 'NUL'); pabs(., 'wi0', bat=0)
    . <- subset(d0('bd0', cache), mtd != 'NUL'); pabs(., 'wd0', bat=0)
    . <- subset(d0('bn0', cache), mtd != 'NUL'); pabs(., 'wn0', bat=0)
    . <- subset(d0('bi0', cache), mtd != 'NUL'); pabs(., 'bi0', bat=1)
    . <- subset(d0('bd0', cache), mtd != 'NUL'); pabs(., 'bd0', bat=1)
    . <- subset(d0('bn0', cache), mtd != 'NUL'); pabs(., 'bn0', bat=1)
}

rpt2 <- function(cache=FALSE)
{
    . <- subset(d0('bz0', cache), mtd %in% c('MNQ', 'GCT', 'BMQ'))
    pabs('n00', FALSE, bias=FALSE, ref=0)
}
