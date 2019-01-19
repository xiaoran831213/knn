library(ggplot2)
library(dplyr)
source('R/hlp.R')

.th <- theme(
    axis.title.x=element_blank(), axis.title.y=element_blank(), 
    strip.text.x = element_text(size=12, face="bold"),
    strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="red", fill="#CCCCFF"),
    legend.title=element_blank(), legend.position='bottom')

## cap the values
.cp <- function(dat, grp, val='val', cap=0.01)
{
    grp <- split(dat, dat[, grp])
    grp <- lapply(grp, function(g)
    {
        v <- g[, val]
        v <- pmin(v, quantile(v, 1-cap/2, na.rm=TRUE))
        v <- pmax(v, quantile(v, 0+cap/2, na.rm=TRUE))
        g[, val] <- v
        g
    })
    dat <- do.call(rbind, grp)
    dat
}

readSIM <- function(sim, cache=FALSE)
{
    rds=paste0(sim, '.rds')
    if(file.exists(rds) && cache)
        agg <- readRDS(rds)
    else
    {
        agg <- lapply(dir(sim, 'rds$', full=TRUE), function(f)
        {
            print(f)
            readRDS(f)
        })
        agg <- do.call(rbind, agg)
        saveRDS(agg, rds)
    }
    invisible(agg)
}

plotErr <- function(sim, out=paste0(sim, '_err.png'))
{
    dat <- readSIM(sim, TRUE)
    dat <- subset(dat, dat=="evl" & key %in% c("ms1", "ms2", "nlk", "cr1", "cr2"))
    dat <- .cp(dat, c("key", "tag"), "val", cap=0.1)

    g <- ggplot(dat, aes(x=mtd, y=val))
    g <- g + geom_boxplot()
    g <- g + facet_grid(key ~ tag, scales="free_y")
    g <- g + .th

    print(out)
    nfx <- length(unique(dat$tag))
    nfy <- length(unique(dat$key))
    ufx <- 19.2 / nfx
    ufy <- ufx / 19.2 * 10.8
    ggsave(out, g, width=19.2, height=min(10.8, ufy * nfy))
    invisible(dat)
}

plotPar <- function(sim, cache=TRUE, out=paste0(sim, '_par.png'))
{
    dat <- readSIM(sim, cache)
    dat <- subset(dat, dat=="par")
    dat <- na.omit(dat)
    dat <- .cp(dat, c("key", "tag"), "val", cap=0.1)

    g <- ggplot(dat, aes(x=key, y=val))
    g <- g + geom_boxplot()
    g <- g + facet_grid(mtd ~ tag, scales="free_y")
    g <- g + .th

    print(out)
    nfx <- length(unique(dat$tag))
    nfy <- length(unique(dat$key))
    ufx <- 19.2 / nfx
    ufy <- 10.8 / 19.2 * ufx
    ggsave(out, g, width=19.2, height=min(10.8, ufy * nfy))
    invisible(dat)
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
        g <- g + .th
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
        g <- g + .th
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
lbl <- function (labels, multi_line = TRUE) 
{
    labels <- label_value(labels, multi_line = multi_line)
    dc <- c(
        `ref`="bold(y)*'='*bold(alpha)+bold(epsilon)",
        `bin`="binomial", `chi`="'chi-square'", `poi`="Poisson", `exp`="exponential",
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
        `n2k`="N*'='*2000", `n4k`="N*'='*4000", `n6k`="N*'='*6000", `n8k`="N*'='*8000",
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

pabs <- function(d, o=NULL, bat=FALSE, xtk=FALSE)
{
    th <- theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(),
                axis.title.y=element_blank(),
                strip.text.x = element_text(size=12, face="bold"),
                strip.text.y = element_text(size=12, face="bold"),
                strip.background = element_rect(colour="red", fill="#CCCCFF"),
                legend.title=element_blank(), legend.position='bottom')
    if(!xtk)
        th <- th + theme(axis.text.x=element_blank())
    
    ## method dictionary
    if(!bat)
        d <- subset(d, !grepl('^B', mtd))

    dc0 <- c(GCT="C0", MNQ="C1", BMQ="C2", BM0="C2", BM1="C2", BM2="C2", BM3="C2", BM4="C2")
    dc2 <- c(GCT="00", MNQ="10", BMQ="30", BM0="20", BM1="21", BM2="22", BM3="23", BM4="24")
    d <- within(d,
    {
        mdc <- dc0[mtd]                 # method category
        mdi <- dc2[mtd]                 # method index
    })
    
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
    g <- ggplot(d, aes(x=mdi, y=val, fill=mdc))
    g <- g + geom_boxplot()
    g <- g + facet_grid(key ~ tag, scales='free_y', labeller=lbl)
    g <- g + .th
    ## g <- g + ggtitle(ttl(d))

    ## fill 3 major categories of methods
    g <- g + scale_fill_manual(
                 breaks=c('C0', 'C1', 'C2'),
                 values=c("#FF7F7F", "#7FFF7F", "#7F7FFF"),
                 labels=c("GCTA   ", "KNN    ", "KNN-Batched"))

    ## x-axis ticks
    g <- g + scale_x_discrete(
                 labels=c("00" = "GCTA",
                          "10" = "KNN",
                          "30" = "B-KNN",
                          "20" = expression(over(KNN, 80)),
                          "21" = expression(over(KNN, 40)),
                          "22" = expression(over(KNN, 20)),
                          "23" = expression(over(KNN, 10))))

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

rpt4 <- function()
{
    a00 <- readSIM('sim/run/a00', FALSE)
    a01 <- readSIM('sim/run/a01', FALSE)
    plotErr('sim/run/a00')
    plotErr('sim/run/a01')
    plotPar('sim/run/a00')
    plotPar('sim/run/a01')
}
rpt0 <- function(cache=TRUE)
{
    . <- subset(d0('bi1', cache), mtd!='NUL' & tag!='ref'); pabs(., 'wi1', bat=0)
    . <- subset(d0('bi1', cache), mtd!='NUL' & tag!='ref'); pabs(., 'bi1', bat=1)
    
    . <- subset(d0('bd1', cache), mtd!='NUL' & tag!='ref'); pabs(., 'wd1', bat=0)
    . <- subset(d0('bd1', cache), mtd!='NUL' & tag!='ref'); pabs(., 'bd1', bat=1)

    . <- subset(d0('bn1', cache), mtd != 'NUL'); pabs(., 'wn1', bat=0)
    . <- subset(d0('bn1', cache), mtd != 'NUL'); pabs(., 'bn1', bat=1)
}

rpt1 <- function(cache=TRUE)
{
    . <- subset(d0('bi0', cache), mtd!='NUL'); pabs(., 'wi0', bat=0)
    . <- subset(d0('bd0', cache), mtd!='NUL'); pabs(., 'wd0', bat=0)
    . <- subset(d0('bn0', cache), mtd!='NUL'); pabs(., 'wn0', bat=0)
    
    . <- subset(d0('bi0', cache), mtd!='NUL'); pabs(., 'bi0', bat=1)
    . <- subset(d0('bd0', cache), mtd!='NUL'); pabs(., 'bd0', bat=1)
    . <- subset(d0('bn0', cache), mtd!='NUL'); pabs(., 'bn0', bat=1)
}

rpt2 <- function(cache=FALSE)
{
    . <- subset(d0('bz0', cache), mtd!='NUL' & tag!='ref'); pabs(., 'bz0', bat=1, xtk=TRUE)
}

## kernel decay test
kdt <- function(cache=TRUE)
{
    d0 <- readSIM('sim/run/kd1', cache)
    d0 <- within(d0, P <- factor(P))
    g <- ggplot(d0, aes(x=P, y=kpa))
    ## g <- g + geom_abline(slope=0, intercept=1, color="red", alpha=.3)
    g <- g + scale_y_log10() 
    g <- g + geom_boxplot()
    g <- g + facet_grid(N ~ K)
    invisible(g)
}
