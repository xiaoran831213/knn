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

    ## dots
    dot <- list(...)
    wtm <- dot$wtm %||% 36              # wall time  (def=36 sec)
    wep <- dot$wep %||% ceiling(nbt)    # wall epoch (def=nbt)
    wit <- dot$wit %||% 1000

    ## message tracks
    tks <- list(msg(~ep, "%04d"), msg(~bt, "%04d"), msg(~mse, "%7.3f"),
                msg(~phi, "%7.3f"), msg(~rtm, "%4.1f"))
    
    ## initial parameters
    itr <- 0L
    rtm <- 0
    ep <- 0
    bt <- 0
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
            . <- sq[seq.int(bt * bsz + 1L, l=bsz) %% N]
            bat <- list(knl=lapply(knl, `[`, ., .), rsp=rsp[.])
        }
        else
        {
            ep <- itr
            bt <- 0
            knl.bat <- knl
            rsp.bat <- rsp
        }
        if(ep > wep) {cat('BMQ: reach max_epoc:', wep, 'e\n'); break}
        t1 <- Sys.time()
        td <- t1 - t0; units(td) <- 'secs'; td <- as.numeric(td)
        
        ## MINQUE on each batch
        bat <- c(bat, with(bat, FUN(rsp, knl, ...)))
        rtm <- rtm + bat$rpt['rtm', 'val'] + td
        par <- bat$par

        ## accumulate an epoch
        eph <- CL(eph, bat)
        if(length(eph) == nbt)
        {
            par <- mat(eph)[, 'nlk']
            hst <- CL(hst, list(ep=ep, rtm=rtm, par=par))
        }
        
        ## record each iteration

        ## message tracks
        if(bt == 0)
        {
            if(ep %% 70 == 0)
                cat(hd.msg(tks), "\n", sep="")
            phi <- par[1]
            mse <- vpd(rsp, knl, par, rt=0, ...)['mse']
            cat(ln.msg(tks), "\n", sep="")
        }

        if(itr > wit) {cat('BMQ: reach max_iter:', wit, 'b\n'); break}
        if(rtm > wtm) {cat('BMQ: reach walltime:', wtm, 's\n'); break}
        itr <- itr + 1L
    }
    ## rearrange history
    tck <- EL1(hst, c('ep', 'rtm'), 'd')
    par <- EL2(hst, 'par')
    hst <- list(tck, par=par)
    
    ## mean solution
    ret <- list(par=mean(par))
    rpt <- with(ret, vpd(rsp, knl, par, ...))
    rpt <- rbind(rtm=DF(key='rtm', val=rtm), rpt)
    
    ## return the history and new parameters
    ret <- c(ret, list(rpt=rpt, hst=hst))
    ret
}
