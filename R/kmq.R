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
kpc.mnq <- function(rsp, knl, rsp.evl=NULL, knl.evl=NULL, bsz=N, ...)
{
    N <- NROW(rsp)                      # sample size
    nbt <- N / bsz                      # number of batches
    
    ## list of history
    hst <- list()

    ## dots
    dot <- list(...)
    wtm <- dot$wtm %||% 36              # wall time  (def=36 sec)
    wep <- dot$wep %||% ceiling(nbt)    # wall epoch (def=nbt)
    max.itr <- dot$max.itr %||% 1000

    ## message tracks
    tks <- list(msg(~ep, "%04d"), msg(~bt, "%04d"), msg(~mse, "%7.3f"),
                msg(~phi, "%7.3f"), msg(~rtm, "%4.1f"))
    
    ## initial parameters
    i <- 0
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
            knl.bat <- lapply(knl, `[`, ix, ix)
            rsp.bat <- rsp[ix]
        }
        else
        {
            ep <- i
            bt <- 0
            knl.bat <- knl
            rsp.bat <- rsp
        }
        if(ep > wep)
        {
            cat('BMQ: reaching max epoch:', wep, '\n');
            break
        }
        t1 <- Sys.time()
        td <- t1 - t0; units(td) <- 'secs'; td <- as.numeric(td)
        
        ## MINQUE on each batch
        bat <- knl.mnq(rsp.bat, knl.bat, cpp=TRUE, ...)
        rtm <- rtm + bat$rpt[1, 2] + td

        ## record each iteration
        hst[[i+1]] <- list(ep=ep, bt=bt, rtm=rtm, par=bat$par, se2=bat$se2)

        ## message tracks
        if(bt == 0)
        {
            par <- bat$par
            phi <- par[1]
            mse <- knl.prd(rsp, knl, par, logged=FALSE, ...)[1, 2]
            cat(ln.msg(tks), "\n", sep="")
        }

        if(i > max.itr) {cat('BMQ: reach max iter:', max.itr, '\n'); break}
        if(rtm > wtm)   {cat('BMQ: reach walltime:', wtm, 'secs\n'); break}
        i <- i + 1
    }
    ## rearrange history
    tck <- EL1(hst, c('ep', 'bt', 'rtm'), 'd')
    par <- EL2(hst, 'par', 'd')
    se2 <- EL2(hst, 'se2', 'd')
    hst <- cbind(tck, par, se2)
    
    ## mean solution
    ret <- list(par=colMeans(par), se2=colMeans(se2))
    rpt <- with(ret, knl.prd(rsp, knl, par, ...))
    rpt <- rbind(rtm=DF(key='rtm', val=rtm), rpt)
    
    ## return the history and new parameters
    ret <- c(ret, list(rpt=rpt, hst=hst))
    ret
}
