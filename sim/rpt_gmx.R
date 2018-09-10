source('sim/rpt.R')

## fetch reports
d0 <- function(name, use.cache=TRUE)
{
    fdr <- file.path('sim/run', name)
    rds <- paste0(fdr, '.rds')
    if(file.exists(rds) && use.cache)
        dat <- readRDS(rds)
    else
    {
        err <- readRpt(fdr, ref=TRUE)
        rtm <- readRtm(fdr, ref=TRUE)
        dat <- within(rbind(rtm, err), {H <- N * R; N <- N * Q})
        dat <- subset(dat, se=-c(Q, R))
        saveRDS(dat, rds)
    }
    dat
}

m0 <- function(name='u20')
{
    img <- file.path('sim/run', name, 'img')
    dat <- d0(name, FALSE)
    for(m in c('mnq'))           # method
    {
        ## errors
        for(a in c('cyh', 'loo', 'mse', 'nlk', 'ssz')) # agg
        {
            . <- c('gct', paste0(m, '.', c('avg', 'whl', 'ssz', a)))
            for(k in c('mse', 'nlk'))   # key
            {
                dt <- subset(dat, mtd %in% . & key == k)
                dt <- within(dt, mtd <- sub('^.*[.]', '', mtd))
                od <- file.path(img, k)

                dir.create(od, FALSE, TRUE)
                fn <- file.path(od, paste0('mbt_', m, '_', a, '.png'))
                print(fn)
                plt(dt, oxy=bsz~sim, ipt=bxp, out=fn)
            }
        }

        ## running time
        {
            . <- c('gct', paste0(m, '.', c('bat', 'whl')))
            k <- 'rtm'
            {
                dt <- subset(dat, mtd %in% . & key == k)
                dt <- within(dt, mtd <- sub('^.*[.]', '', mtd))
                dt <- within(dt, mtd <- relevel(as.factor(mtd), 'gct'))
                od <- file.path(img, k)

                dir.create(od, FALSE, TRUE)
                fn <- file.path(od, paste0('mbt_', m, '.png'))
                print(fn)
                plt(dt, oxy=bsz~sim, ipt=bxp, out=fn)
            }
        }
    }
}

## whole sample single cohort, GCTA vs. MINQUE
m1 <- function(name='u13')
{
    img <- file.path('sim/run', name, 'img')
    d <- d0(name, FALSE)
    
    ## MSE, NLK, and time, hide svc
    d <- subset(d, key %in% c('mse', 'cyh', 'rtm'), -svc)

    ## 1st order data, and 1st order polynomial kernel
    dir.create(img, FALSE, TRUE)

    ## 1st order data, 2nd order polynomial kernel
    plt(subset(d, sim == 'I(p1)~p2, x~a1'), oxy=.~key, ipt=bxp, out=file.path(img, 'ukb_whl_p01.png'))

    ## 2st order data, 2nd order polynomial kernel
    plt(subset(d, sim == 'I(p1)~p2, x~a2'), oxy=.~key, ipt=bxp, out=file.path(img, 'ukb_whl_p02.png'))

    ## plt(d, oxy=sim~key, ipt=bxp, out=file.path(img, 'whl.png'))
}

## batched traning, single cohort test
m2 <- function(name='u20')
{
    img <- file.path('sim/run', name, 'img')
    d <- d0(name, FALSE)
    dir.create(img, FALSE, TRUE)

    ## GCTA vs. MINQUE, whole sample and naive batch
    d <- subset(d, mtd == 'gct' | grepl('^mnq[.](whl|ssz|bat)$', mtd))
    d <- within(d, mtd <- sub('[.]ssz$', '.bat', mtd))

    ## whole sample means batch size == sample size
    d <- within(d, mtd[mtd == 'mnq.whl'] <- 'whl')
    d <- within(d, mtd[mtd == 'mnq.bat'] <- bsz[mtd == 'mnq.bat'])
    d <- within(d, mtd[mtd == 'mnq.ssz'] <- 'bat')

    ## MSE, NLK, and time, hide svc
    d <- subset(d, key %in% c('mse', 'nlk'), -svc)

    ## 1st order data, and 2st order polynomial kernel
    plt(subset(d, sim == 'I(p1)~p2, x~a1'), oxy=.~key, ipt=bxp, out=file.path(img, 'ukb_bat_p01.png'))
                        
    ## 2nd order data, and 2nd order polynomial kernel
    plt(subset(d, sim == 'I(p1)~p2, x~a2'), oxy=.~key, ipt=bxp, out=file.path(img, 'ukb_bat_p02.png'))
}

## whole sample, timing by sample size
t1 <- function(name='t10')
{
    library(ggplot2)
    library(dplyr)

    ## output path
    img <- file.path('sim/run', name, 'img')
    
    ## report files
    fns <- dir(file.path('sim/run', name), 'rds$', ful=TRUE)
    rpt <- list()
    for(f in fns)
    {
        print(f)
        rpt <- c(rpt, list(readRDS(f)))
    }
    rpt <- do.call(rbind, rpt)

    ## timing related only
    dat <- as_tibble(rpt)
    dat <- dat %>% filter(dat=='dvp', key=='rtm', N<=8000) %>% sel(-dat, -key, -seed)
    dat <- dat %>% filter(yks=='p2')
    grp <- dat %>% group_by_at(vars(N:mtd))
    agg <- grp %>% summarize(val=mean(val), std=sd(val)) %>% ungroup

    ## plot
    g <- ggplot(agg, aes(x=N, y=val))
    g <- g + geom_line(aes(color=mtd))
    ## g <- g + facet_wrap(~yks)

    dir.create(img, FALSE, TRUE)
    ggsave(file.path(img, 'ukb_rtm.png'), width=5, height=4)
    g
}


## quick explore of NLK
prob <- function(name='u12')
{
    library(ggplot2)
    library(dplyr)
    d <- as_tibble(d0(name, FALSE))
    d <- filter(d, key %in% c('cyh', 'mse', 'rtm'))

    g <- ggplot(d, aes(x=mtd, y=val))
    g <- g + geom_boxplot()
    g <- g + ylim(0, 1)
    g <- g + facet_grid(sim~key)
    g <- g + ggtitle(ttl(d))

    ## write to PNG
    f <- paste0(file.path('~', name), '.png')
    dir.create(basename(f), FALSE, TRUE)
    print(f)
    ggsave(f, g, width=10, height=10)

    ## return the plot
    g
}
