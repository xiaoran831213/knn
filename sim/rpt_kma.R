source('sim/rpt.R')

m3 <- function()
{
    rds <- 'sim/run/s3x.rds'
    if(file.exists(rds))
        da <- readRDS(rds)
    else
    {
        s0 <- subset(readRpt('sim/run/s30'))
        s1 <- subset(readRpt('sim/run/s31'))
        s2 <- subset(readRpt('sim/run/s32'))
        s3 <- subset(readRpt('sim/run/s33'))
        da <- rbind(s0, s1, s2, s3)
        saveRDS(da, rds)
    }

    ## s0 <- within(subset(s0, grepl('[.]mnq$', mtd)), {use <- 'mnq'; mtd <- sub('[.]mnq$', '', mtd)})
    ## s1 <- within(subset(s1, grepl('[.]mnq$', mtd)), {use <- 'mnq'; mtd <- sub('[.]mnq$', '', mtd)})
    ## s2 <- within(subset(s2, grepl('[.]mnq$', mtd)), {use <- 'mnq'; mtd <- sub('[.]mnq$', '', mtd)})
    ## s3 <- within(subset(s3, grepl('[.]mnq$', mtd)), {use <- 'mnq'; mtd <- sub('[.]mnq$', '', mtd)})
    ## plt(s0, oxy=.~het, ipt=bxp, out='rpt/2018_07_08/km3_mnq_s00.png')
    ## plt(s1, oxy=.~het, ipt=bxp, out='rpt/2018_07_08/km3_mnq_s01.png')
    ## plt(s2, oxy=.~het, ipt=bxp, out='rpt/2018_07_08/km3_mnq_s02.png')
    ## plt(s3, oxy=.~het, ipt=bxp, out='rpt/2018_07_08/km3_mnq_s03.png')

    da <- subset(da, key!='rtm')
    d1 <- within(subset(da, grepl('[.]mnq$', mtd)), {use <- 'mnq'; mtd <- sub('[.]mnq$', '', mtd)})
    d2 <- within(subset(da, grepl('[.]rop$', mtd)), {use <- 'rop'; mtd <- sub('[.]rop$', '', mtd)})
    plt(d1, oxy=N~het, ipt=bxp, axi=val~mtd, out='rpt/2018_07_08/km3_mnq.png')
    plt(d2, oxy=N~het, ipt=bxp, axi=val~mtd, out='rpt/2018_07_08/km3_rop.png')

    ## cohort size is 250
    db <- subset(da, N==250 & grepl('^st', sim) & het %in% c(0.0, 0.50, 1.0))
    b1 <- within(subset(db, grepl('[.]mnq$', mtd)), {use <- 'mnq'; mtd <- sub('[.]mnq$', '', mtd)})
    b2 <- within(subset(db, grepl('[.]rop$', mtd)), {use <- 'rop'; mtd <- sub('[.]rop$', '', mtd)})
    b1 <- subset(b1, mtd %in% c('avg', 'ssz', 'whl'))
    b2 <- subset(b2, mtd %in% c('avg', 'ssz', 'whl'))

    ## name change for methods
    dict <- c(avg='avg', ssz='meta', whl='mega')
    b1 <- within(b1, mtd <- dict[mtd])
    b2 <- within(b2, mtd <- dict[mtd])

    plt(b1, oxy=key~het, ipt=bxp, axi=val~mtd, out='rpt/kma/img/jsm_mm3_mnq.png')
    plt(b2, oxy=key~het, ipt=bxp, axi=val~mtd, out='rpt/kma/img/jsm_mm3_rop.png')
}
