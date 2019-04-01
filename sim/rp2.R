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
.cp <- function(dat, grp, val='val', cap=0.01, mtd=c('both', 'upper', 'lower'))
{
    grp <- split(dat, dat[, grp])
    mtd <- match.arg(mtd, c('both', 'upper', 'lower'))
    grp <- lapply(grp, function(g)
    {
        v <- g[, val]
        if(mtd == 'upper')
            v <- pmin(v, quantile(v, 1-cap, na.rm=TRUE))
        else if(mtd == 'lower')
            v <- pmax(v, quantile(v, 0+cap, na.rm=TRUE))
        else
        {
            v <- pmin(v, quantile(v, 1-cap/2, na.rm=TRUE))
            v <- pmax(v, quantile(v, 0+cap/2, na.rm=TRUE))
        }
        g[, val] <- v
        g
    })
    dat <- do.call(rbind, grp)
    dat
}

readSIM <- function(sim, cache=TRUE)
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

plotErr <- function(sim, out=paste0(sim, '_err.png'), cap=0.05, ...)
{
    dat <- readSIM(sim, TRUE)
    dat <- subset(dat,
                  dat=="evl" & key %in% c("yel", "ycl", "nlk"))
                  ## dat=="evl" & key %in% c("yeh", "ych")|                  
                  ## dat=="dvp" & key %in% c("hsq"))
    dat$key <- as.character(dat$key)
    grp <- dat[, c("seed", "key", "tag", "dat")]
    
    ## minus NUL effect
    dat <- by(dat, grp, function(g)
    {
        m <- g$mtd != "NUL"
        nul <- g[!m, "val"]
        g <- g[m, ]
        if(!g$key[1] %in% c("ycl", "ych", "hsq"))
            g$val <- g$val / nul
        g
    }, simplify=FALSE)
    dat <- do.call(rbind, dat)
    dat <- .cp(dat, c("key", "tag", "dat"), "val", cap=0.05, ...)
    dat <- within(dat, val <- pmin(1.1, val))
    ## dat <- subset(dat, !tag %in% c("bas", "ref"))

    g <- ggplot(dat, aes(x=mtd, y=val))
    g <- g + geom_boxplot()
    g <- g + facet_grid(key ~ tag, scales="free_y", labeller=lbl)
    g <- g + .th
    
    nfy <- length(unique(dat$key))
    ufy <- 10 / nfy
    nfx <- length(unique(dat$tag))
    ufx <- 19 / nfx
    if(ufx / ufy < 19 / 10)
        ufy <- ufx / 19 * 10
    else
        ufx <- ufy / 10 * 19
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy))
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


## labeller
lbl <- function (labels, multi_line = TRUE) 
{
    labels <- label_value(labels, multi_line = multi_line)
    dc <- c(
        `ref`="bold(y)*'='*bold(g)+bold(epsilon)",
        `bin`="binomial", `chi`="'chi-square'", `poi`="Poisson", `exp`="exponential",
        `hy1`="bold(y)*'='*over('|'*bold(g)*'|', 1 + '|'*bold(g)*'|') + bold(epsilon)",
        `rc1`="bold(y)*'='*bold('|'*g*'|')*e^bold('|'*g*'|') + bold(epsilon)",
        `g^2`="bold(y)*'='*bold(g)^2 + bold(epsilon)",
        `g^3`="bold(y)*'='*bold(g)^3 + bold(epsilon)",
        `g:2`="'2-way interaction'",
        `g:3`="'3-way interaction'",
        `bs0`="hat(symbol(sigma))[0]^2 - symbol(sigma)[0]^2",
        `nlk`="NLK(bold(y), hat(bold(theta)))",
        `yel`="MSE(bold(y), hat(bold(y)))",
        `ycl`="COR(bold(y), hat(bold(y)))",
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

main.plot <- function()
{
    plotErr('sim/run/a09')
    plotErr('sim/run/n11')
    ## plotErr('sim/run/sm0')
    plotErr('sim/run/bn0')
    invisible(NULL)
}
