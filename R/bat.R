## General Batched Trainer

batch.mask <- function(rsp, nbt=5, bpt=0, ...)
{
    N <- NROW(rsp)
    msk <- replicate(ceiling(N / nbt), sample.int(nbt))[seq.int(N)]
    ## cnt <- sapply(seq(nbt), function(x) sum(msk == x))
    ## names(cnt) <- seq(nbt)
    ## print(list(cnt=cnt))
    if(bpt==1)
    {
        seq <- sort(rsp, index.return=TRUE)$ix
        qes <- sort(seq, index.return=TRUE)$ix
        msk <- msk[qes]
    }
    else
        msk <- sample(msk)

    ## print(list(var=var(rsp)))
    ## for(i in seq.int(nbt))
    ## {
    ##     y <- rsp[msk == i]
    ##     print(var(y))
    ## }
    invisible(msk)
}

#' @title General Batched Trainer
#' @param fun function the non-batched algorithm
#' @param nbt number of batches (def=5)
#' @param wtm numeric (in ...) walltime in hours
#' @param wep integer (in ...) wall epoch count
#'
#' @return a list of training history, the averaged parameters
GBT <- function(FUN, rsp, knl, xmx=NULL, bsz=100, ...)
{
    N <- NROW(rsp)                      # sample size
    nbt <- as.integer(N / bsz)
    ## list of history
    hst <- list()

    ## dots
    dot <- list(...)
    wtm <- dot$wtm %||% 36              # wall time  (def=36 sec)
    wep <- dot$wep %||% 1               # wall epoch (def=1)
    wit <- dot$wit %||% 1000
    pss <- dot$pss %||% 0               # pass on estimates?
    ini <- dot$ini
    bpt <- dot$bpt %||% 0

    ## message tracks
    tks <- list(msg(~ep, "%04d"),
                msg(~mse, "%7.3f"), msg(~nlk, "%7.3f"),
                msg(~eps, "%7.3f"), msg(~rtm, "%4.1f"))
    
    ## initial parameters
    itr <- 1; rtm <- 0; ep <- 1; bt <- 1
    while(TRUE)
    {
        ## create batches
        t0 <- Sys.time()
        if(bt == 1L)                    # new batch
        {
            bmk <- batch.mask(rsp, nbt, bpt) # batch mask
            eph <- list()
        }
        
        ## get batch data
        . <- bmk == bt
        bat <- list(knl=lapply(knl, `[`, ., .), rsp=rsp[.])
        bat$xmx <- if(is.null(xmx)) NULL else xmx[., , drop=FALSE]
        

        ## get batch estimates
        bat <- c(bat, with(bat, FUN(rsp, knl, xmx, ini, ...)))
        ini <- if(pss) bat$par else NULL # pass on initial values?

        t1 <- Sys.time()
        td <- t1 - t0; units(td) <- 'secs'; td <- as.numeric(td)
        rtm <- rtm + td

        ## gather information
        msg <- c(ep=ep, bt=bt, rtm=rtm, bat[c('par', 'se2')], bsz=bat$rpt['ssz', 'val'])
        
        ## accumulate an epoch, calculate training statistics
        eph <- CL(eph, bat[c('par', 'se2')])
        bt <- bt + 1
        if(bt == nbt)                   # end of an epoch?
        {
            par  <- if(pss) bat[['par']] else par <- mean(eph %$% 'par')
            eps <- par['EPS']
            rpt <- vpd(rsp, knl, xmx, par, rt=0, ...)
            mse <- rpt['yel']
            nlk <- rpt['nlk']
            msg <- c(msg, rpt)          # additional info
            ## append message to STDOUT
            cat(ln.msg(tks), "\n", sep="")
            ep <- ep + 1L               # next epoch
            bt <- 1L                    # next epoch
        }

        ## append message to history
        hst <- CL(hst, msg)
        
        if(rtm > wtm) {cat('BMQ: reach walltime:', wtm, 's\n'); break}
        if(ep > wep) {cat('BMQ: reach max_epoc:', wep, 'e\n'); break}
    }
    ## mean solution
    ## par <- mean(hst %$% 'par')
    par <- do.call(rbind, hst %$% 'par')
    ## par <- apply(par, 2L, geom.mean)
    par <- apply(par, 2L, mean)

    ## assessment
    t0 <- Sys.time(); rpt <- vpd(y=rsp, K=knl, W=par); t1 <- Sys.time()
    td <- t1 - t0; units(td) <- 'secs'; td <- as.numeric(td)
    rtm <- rtm + td
    rpt <- rbind(rpt, rtm=DF(key='rtm', val=rtm))

    ## return
    ret <- list(par=par, rpt=rpt) #, hst=hst)
    ret
}
