source('sim/rpt.R')
## create plot for JSM
m1.jsm <- function()
{
    rds <- 'sim/run/jsm_bt1.rds'
    if(file.exists(rds))
        d <- readRDS(rds)
    else
    {
        f <- c('sim/run/b10', 'sim/run/b11', 'sim/run/b12', 'sim/run/b13')
        rtm <- readRtm(f)
        err <- readRpt(f)
        d <- rbind(rtm, err)
        saveRDS(d, rds)
    }
    d <- subset(d, het==0, -c(het, bsz, wep))

    gct <- subset(d, sim=='I(ga)~p1' & mtd=='gct')
    gct$mtd <- 'GCTA'
    mq1 <- subset(d, sim=='I(ga)~p1' & grepl('(mnq|nlk).whl', mtd))
    mq1$mtd <- 'MNQ1'
    mq2 <- subset(d, sim=='I(ga)~p2' & grepl('(mnq|nlk).whl', mtd))
    mq2$mtd <- 'MNQ2'
    rpt <- rbind(gct, mq1, mq2)
    plt(rpt, oxy=.~key, ipt=bxp) ##, out='rpt/kma/img/vcm_bmk_gau.png')
}

m2.jsm <- function()
{
    rds <- 'sim/run/jsm_bt1.rds'
    if(file.exists(rds))
        d <- readRDS(rds)
    else
    {
        f <- c('sim/run/b10', 'sim/run/b11', 'sim/run/b12', 'sim/run/b13')
        rtm <- readRtm(f)
        err <- readRpt(f)
        d <- rbind(rtm, err)
        saveRDS(d, rds)
    }
    d <- subset(d, het==0, -c(het, bsz, wep))
    gct <- subset(d, sim=='I(p2)~p1' & mtd=='gct')
    gct$mtd <- 'GCTA'
    mq1 <- subset(d, sim=='I(p2)~p1' & grepl('(mnq|nlk).whl', mtd))
    mq1$mtd <- 'MNQ1'
    mq2 <- subset(d, sim=='I(p2)~p2' & grepl('(mnq|nlk).whl', mtd))
    mq2$mtd <- 'MNQ2'
    rpt <- rbind(gct, mq1, mq2)
    plt(rpt, oxy=.~key, ipt=bxp) ## , out='rpt/kma/img/vcm_bk2_gau.png')
}


m0.jsm <- function()
{
    rds <- 'sim/run/jsm_mt0.rds'
    if(file.exists(rds))
    {
        d <- readRDS(rds)
    }
    else
    {
        f <- c('sim/run/s00', 'sim/run/s01')
        rtm <- readRtm(f)
        err <- readRpt(f)
        d <- rbind(rtm, err)
        saveRDS(d, rds)
    }

    ## dc <- rep(c(avg='avg', nlk='validity', loo='validity', ssz='precision', whl='mega'), e=2)
    ## names(dc) <- paste(c('mnq', 'mle'), names(dc), sep='.')
    r <- d %>% group_by_at(vars(N:key, sim)) %>% summarize(val=mean(val)) %>% ungroup
    r <- r %>% filter(key=='nlk') %>% select(het, mtd, val, sim) %>%
        mutate(agg=sub('^.*[.]', '', mtd), alg=sub('[.].*$', '', mtd)) %>%
        select(-mtd)
    r <- r %>% group_by(alg, het) %>% mutate(ref=val/nth(val, 6)) %>% filter(ref<1.0)

    ## r %>% filter(sim=='I(p2+ga)~p2+ga') %>% as.data.frame
    d <- subset(d, sim=='I(p2+ga)~p2+ga' & key != 'rtm')
    d <- subset(d, het %in% c(0.0, 0.5, 1.0))

    ## d <- within(d, mtd <- dc[mtd])
    for(m in c('mle', 'mnq'))           # method
    {
        for(a in c('cyh', 'loo', 'mse', 'nlk'))
        {
            . <- paste0(m, '.', c('avg', 'whl', 'ssz', a))
            print(.)
            ## dev.new()
            dt <- within(subset(d, mtd %in% .), mtd <- sub('^.*[.]', '', mtd))
            
            fn <- paste0('rpt/2018_08_09/img/met_', m, '_', a, '_mse.png')
            plt(subset(dt, key=='mse'), oxy=.~het, ipt=bxp, out=fn)
            fn <- paste0('rpt/2018_08_09/img/met_', m, '_', a, '_nlk.png')
            plt(subset(dt, key=='nlk'), oxy=.~het, ipt=bxp, out=fn)
        }
    }
    ## plt(d, oxy=het~key, ipt=bxp, out='rpt/kma/img/met_het_mnq_nlk.png')
    r
}

## single cohort
c0.jsm <- function()
{
    rds <- 'sim/run/c0x.rds'
    if(file.exists(rds))
    {
        print('read from rds')
        d <- readRDS(rds)
    }
    else
    {
        f <- c('sim/run/c00')
        rtm <- readRtm(f)
        err <- readRpt(f)
        d <- rbind(rtm, err)
        saveRDS(d, rds)
    }

    ## PLY
    gct <- subset(d, sim=='I(p2)~p1' & mtd=='gct')
    gct$mtd <- 'GCTA'

    mq1 <- subset(d, sim=='I(p2)~p1' & mtd=='mnq.whl'); mq1$mtd <- 'MNQ1'
    mq2 <- subset(d, sim=='I(p2)~p2' & mtd=='mnq.whl'); mq2$mtd <- 'MNQ2'
    ml1 <- subset(d, sim=='I(p2)~p1' & mtd=='mle.whl'); ml1$mtd <- 'MLE1'
    ml2 <- subset(d, sim=='I(p2)~p2' & mtd=='mle.whl'); ml2$mtd <- 'MLE2'
    rpt <- rbind(gct, mq1, mq2, ml1, ml2)

    plt(rpt, oxy=.~key, ipt=bxp, out='rpt/kma/img/vcm_pl2_all.png')
    plt(subset(rpt, grepl('^(MLE|GCT)', mtd)), oxy=.~key, ipt=bxp, out='rpt/kma/img/vcm_ply_mle.png')
    plt(subset(rpt, grepl('^(MNQ|GCT)', mtd)), oxy=.~key, ipt=bxp, out='rpt/kma/img/vcm_ply_mnq.png')

    ## GAU
    gct <- subset(d, sim=='I(ga)~p1' & mtd=='gct')
    gct$mtd <- 'GCTA'

    mq1 <- subset(d, sim=='I(ga)~p1' & mtd=='mnq.whl'); mq1$mtd <- 'MNQ1'
    mq2 <- subset(d, sim=='I(ga)~p2' & mtd=='mnq.whl'); mq2$mtd <- 'MNQ2'
    ml1 <- subset(d, sim=='I(ga)~p1' & mtd=='mle.whl'); ml1$mtd <- 'MLE1'
    ml2 <- subset(d, sim=='I(ga)~p2' & mtd=='mle.whl'); ml2$mtd <- 'MLE2'
    rpt <- rbind(gct, mq1, mq2, ml1, ml2)

    plt(rpt, oxy=.~key, ipt=bxp, out='rpt/kma/img/vcm_gau_all.png')
    plt(subset(rpt, grepl('^(MLE|GCT)', mtd)), oxy=.~key, ipt=bxp, out='rpt/kma/img/vcm_gau_mle.png')
    plt(subset(rpt, grepl('^(MNQ|GCT)', mtd)), oxy=.~key, ipt=bxp, out='rpt/kma/img/vcm_gau_mnq.png')
}
