## Kernel MINQUE (KMQ)

#' @title gradient descent for kernel deep net
#' @param par vector initial values of parameters;
#' @param ctx list context of the training task, such as input features and lables;
#' @param usz integer number of sampling to calculate the losses and derivatives;
#' @param bsz integer batch size of stochastic gradient descent;
#' @param mmt numeric momment of stochastic gradient descent, between 0 and 1;
#' @param lr0 numeric the initial learning rate;
#' @param wtm numeric (in ...) walltime in hours
#' @param esz integer (in ...) number of samples to draw in for error estimation;
#' also collect training and evaluation statistics such as L1 and L2 error.
kpc.mnq <- function(rsp.dvp, knl.dvp, rsp.evl=NULL, knl.evl=NULL, bsz=N, ...)
{
    N <- NROW(knl.dvp[[1]])             # sample size
    nbt <- N / bsz                      # number of batches
    
    ## list of history
    hst <- list()

    ## dots
    dot <- list(...)
    wtm <- dot$wtm %||% 36              # wall time  (def=36 sec)
    wep <- dot$wep %||% ceiling(nbt)    # wall epoch (def=nbt)
    max.itr <- dot$max.itr %||% 1000

    ## header of tracks
    hdr <- c(
        eph=sprintf("%4s", 'eph'),         # epoch
        bat=sprintf("%4s", 'bat'),         # batch
        mse.bat=sprintf("%7s", 'mse.bat'), # mse.bat
        mse.dvp=sprintf("%7s", 'mse.dvp'), # mse.dvp
        mse.evl=sprintf("%7s", 'mse.evl'), # mse.evl
        phi=sprintf('%7s', 'phi'),         # phi
        rtm=sprintf('%4s', 'rtm'))         # seconds used
    hdr <- paste(hdr, collapse=' ')

    ## initial parameters
    i <- 0
    sq <- seq.int(N)
    hst.mse <- list()
    hst.par <- list()
    hst.se2 <- list()
    hst.num <- list()
    rtm <- 0
    while(TRUE)
    {
        ## create batches
        t0 <- Sys.time()
        if(nbt > 1)
        {
            ep <- as.integer((i * bsz) / N) + 1 # epoch count
            bt <- (i * bsz) %% N / bsz          # batch count
            if(bt == 0)                         # 1st batch?
                sq <- sample.int(N)             # permutation

            ix <- sq[seq.int(bt * bsz, l=bsz) %% N]
            knl.bat <- lapply(knl.dvp, `[`, ix, ix)
            rsp.bat <- rsp.dvp[ix]
        }
        else
        {
            ep <- i
            bt <- 0
            knl.bat <- knl.dvp
            rsp.bat <- rsp.dvp
        }
        t1 <- Sys.time()
        td <- t1 - t0; units(td) <- 'secs'; td <- as.numeric(td)
        
        ## MINQUE on each batch
        bat.ret <- knl.mnq(rsp.bat, knl.bat, cpp=TRUE, ...)
        rtm <- rtm + bat.ret$rpt[1, 2] # + td
        par <- bat.ret$par
        se2 <- bat.ret$se2
        phi <- par[1]
        bat.rpt <- bat.ret$rpt[-1, ]    # keep errors only
        dvp.rpt <- knl.mnq.evl(rsp.dvp, knl.dvp, par, ...)
        evl.rpt <- knl.mnq.evl(rsp.evl, knl.evl, par, ...)
        mse <- list(bat=bat.rpt[1, 2], dvp=dvp.rpt[1, 2], evl=evl.rpt[1, 2])

        ## record each iteration
        hst.num[[i+1]] <- list(i=i, ep=ep, bt=bt, rtm=rtm)
        hst.par[[i+1]] <- par
        hst.se2[[i+1]] <- se2
        hst.mse[[i+1]] <- mse

        ## update learning rate
        if(i %% 80 == 0)
            cat(hdr, '\n', sep='')
        
        ## message tracks
        msg <- c(
            eph=sprintf("%04X", as.integer(ep)),   # epoch
            bat=sprintf("%04X", as.integer(bt)),   # batch
            mse.bat=sprintf("%7.3f", mse$bat),     # mse.bat
            mse.dvp=sprintf("%7.3f", mse$dvp),     # mse.dvp
            mse.evl=sprintf("%7.3f", mse$evl),     # mse.evl
            phi=sprintf('%7.3f', phi),             # phi
            rtm=sprintf('%4.1f', rtm))             # seconds used
        cat(paste(msg, collapse=' '), '\n', sep='')
        
        if(i > max.itr)
        {
            cat('BMQ: reaching max iter:', max.itr, '\n')
            break
        }
        if(ep > wep)
        {
            cat('BMQ: reaching max iter:', wep, '\n')
            break
        }
        if(rtm > wtm)
        {
            cat('BMQ: reaching walltime:', wtm, 'hour(s)\n')
            break
        }
        i <- i + 1
    }

    ## rearrange history
    hst.num <- do.call(rbind, lapply(hst.num, unlist))
    hst.mse <- do.call(rbind, lapply(hst.mse, unlist))
    hst.par <- do.call(rbind, lapply(hst.par, unlist))
    hst.se2 <- do.call(rbind, lapply(hst.se2, unlist))
    hst <- DF(hst.num, mse=hst.mse, phi=hst.par[, 1])
    
    ## mean parameter solution
    par <- apply(hst.par, 2, mean)
    se2 <- apply(hst.se2, 2, mean)
    rpt <- knl.mnq.evl(rsp.dvp, knl.dvp, par, ...)
    rpt <- rbind(DF(key='rtm', val=tail(hst$rtm, 1)), rpt)

    ## return the history and new parameters
    ret <- list(rpt=rpt, par=par, se2=se2, hst=hst)
    ret
}