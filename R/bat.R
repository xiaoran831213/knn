## General Batched Trainer

#' @title General Batched Trainer
#' @param fun function the non-batched algorithm
#' @param bsz integer batch size
#' @param wtm numeric (in ...) walltime in hours
#' @param wep integer (in ...) wall epoch count
#'
#' @return a list of training history, the averaged parameters
GBT <- function(FUN, rsp, knl, bsz=NROW(rsp), ...)
{
    N <- NROW(rsp)                      # sample size
    nbt <- N / bsz                      # number of batches
    
    ## list of history
    hst <- list()
    acc <- list()

    ## dots
    dot <- list(...)
    wtm <- dot$wtm %||% 36              # wall time  (def=36 sec)
    wep <- dot$wep %||% 1               # wall epoch (def=1)
    wit <- dot$wit %||% 1000

    ## message tracks
    tks <- list(msg(~ep, "%04d"),
                msg(~mse, "%7.3f"), msg(~nlk, "%7.3f"), msg(~loo, "%7.3f"),
                msg(~phi, "%7.3f"), msg(~rtm, "%4.1f"))
    
    ## initial parameters
    itr <- 0; rtm <- 0; ep <- 0; bt <- 0
    while(TRUE)
    {
        ## create batches
        t0 <- Sys.time()
        if(nbt > 1)
        {
            ep <- as.integer((itr * bsz) / N) + 1 # epoch count
            bt <- (itr * bsz) %% N / bsz          # batch count
            if(bt == 0)                           # 1st batch?
            {
                sq <- sample.int(N)               # permutation
                eph <- list()
            }
            . <- sq[seq.int(bt * bsz + 1L, l=bsz) %% (N + 1)]
            bat <- list(knl=lapply(knl, `[`, ., .), rsp=rsp[.])
        }
        else
        {
            ep <- itr + 1
            bt <- 0
            knl.bat <- knl
            rsp.bat <- rsp
        }
        bt <- bt + 1
        if(ep > wep) {cat('BMQ: reach max_epoc:', wep, 'e\n'); break}
        t1 <- Sys.time()
        td <- t1 - t0; units(td) <- 'secs'; td <- as.numeric(td)
        
        ## MINQUE on each batch
        ssz <- nrow(bat)
        bat <- c(bat, with(bat, FUN(rsp, knl, ...)))
        rtm <- rtm + bat$rpt['rtm', 'val'] + td

        ## gather information
        msg <- c(ep=ep, bt=bt, rtm=rtm, bat[c('par', 'se2')], bsz=bat$rpt['ssz', 'val'])
        
        ## accumulate an epoch, calculate training statistics
        eph <- CL(eph, bat[c('par', 'se2')])
        if(length(eph) == nbt)
        {
            par <- mean(eph %$% 'par')
            phi <- par[1]
            rpt <- vpd(rsp, knl, par, rt=0, ...)
            mse <- rpt['mse']
            nlk <- rpt['nlk']
            loo <- rpt['loo']
            msg <- c(msg, rpt)
            ## append message to STDOUT
            cat(ln.msg(tks), "\n", sep="")
        }

        ## append message to the history
        hst <- CL(hst, msg)
        
        if(itr > wit) {cat('BMQ: reach max_iter:', wit, 'b\n'); break}
        if(rtm > wtm) {cat('BMQ: reach walltime:', wtm, 's\n'); break}
        itr <- itr + 1L
    }
    
    ## mean solution
    par <- mean(hst %$% 'par')

    ## assessment
    t0 <- Sys.time(); rpt <- vpd(y=rsp, K=knl, W=par); t1 <- Sys.time()
    td <- t1 - t0; units(td) <- 'secs'; td <- as.numeric(td)
    rtm <- rtm + td
    rpt <- rbind(rpt, rtm=DF(key='rtm', val=rtm))

    ## return
    ret <- list(par=par, rpt=rpt) #, hst=hst)
    ret
}
