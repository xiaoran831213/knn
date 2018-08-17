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
    wep <- dot$wep %||% ceiling(nbt)    # wall epoch (def=nbt)
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
            . <- sq[seq.int(bt * bsz + 1L, l=bsz) %% N]
            bat <- list(knl=lapply(knl, `[`, ., .), rsp=rsp[.])
        }
        else
        {
            ep <- itr + 1
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

        ## processing at the end of an epoch
        eph <- CL(eph, bat)
        if(length(eph) == nbt)
        {
            par <- mat(eph)
            phi <- par[1, 'ssz']
            rpt <- vpd(rsp, knl, par, rt=0, ...)
            mse <- rpt['mse']
            nlk <- rpt['nlk']
            loo <- rpt['loo']
            cat(ln.msg(tks), "\n", sep="")
            hst <- CL(hst, c(ep=ep, rtm=rtm, rpt, list(par=par)))
        }
        acc <- CL(acc, bat)
        
        ## append message track on each new epoch
        if(bt == 0 && ep %% 40 == 1)
            cat(hd.msg(tks), "\n", sep="")

        if(itr > wit) {cat('BMQ: reach max_iter:', wit, 'b\n'); break}
        if(rtm > wtm) {cat('BMQ: reach walltime:', wtm, 's\n'); break}
        itr <- itr + 1L
    }
    ## rearrange history
    tck <- EL1(hst, c('ep', 'rtm', 'mse', 'loo', 'nlk', 'cyh'), 'd')

    
    ## mean solution
    ## par <- EL2(hst, 'par')
    ## par <- mean(par)
    par <- mat(acc)
    par <- sapply(colnames(par), function(i) par[, i], simplify=FALSE)
    rpt <- lapply(par, vpd, y=rsp, K=knl)

    ## return
    ret <- list(par=par, rpt=rpt, rtm=rtm, hst=tck)
    ret
}
