## Kernel MINQUE (KMQ)

## Training Track
msg <- function(obj, fmt=NULL, hdr=NULL, ...)
{
    if(is.function(obj))
    {
        ret <- function()
        {
            sprintf(fmt, obj(...))
        }
    }
    else if(inherits(obj, 'formula'))
    {
        ret <- function()
        {
            e <- environment(obj)
            sprintf(fmt, eval(attr(terms(obj), 'variables'), e)[[1]])
        }
    }
    else
    {
        ret <- function()
        {
            sprintf(fmt, obj)
        }
    }

    if(is.null(fmt))
        fmt <- "%8s"

    if(is.null(hdr))
        hdr <- deparse(substitute(obj))
    hdr <- sub("^~", "", hdr)
    
    len <- as.integer(sub("^%[+-]?([0-9]*)[.]?([0-9])?[A-z]$", "\\1", fmt))
    hdr <- sprintf(paste0("% ", len, "s"), hdr)
    structure(ret, class=c('msg', 'function'), hdr=hdr)
}

## is
is.msg <- function(.) 'msg' %in% class(.)

## Header of the tracks
hd.msg <- function(...)
{
    d <- Filter(is.msg, unlist(list(...)))
    d <- sapply(d, function(.)
    {
        h <- attr(., "hdr")
    })
    paste(d, collapse=" ")
}

## A line of the tracks
ln.msg <- function(...)
{
    d <- Filter(is.msg, unlist(list(...)))
    d <- sapply(d, do.call, list())
    paste(d, collapse=" ")
}


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

    ## message tracks
    tks <- list(msg(~ep, "%04d"), msg(~bt, "%04d"), msg(~mse, "%7.3f"),
                msg(~phi, "%7.3f"), msg(~rtm, "%4.1f"))
    
    ## initial parameters
    i <- 0
    hst.mse <- list()
    hst.se2 <- list()
    hst.par <- list()
    hst.num <- list()
    hst.err <- list()
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
        if(ep > wep)
        {
            cat('BMQ: reaching max epoch:', wep, '\n');
            break
        }
        t1 <- Sys.time()
        td <- t1 - t0; units(td) <- 'secs'; td <- as.numeric(td)
        
        ## MINQUE on each batch
        bat.ret <- knl.mnq(rsp.bat, knl.bat, cpp=TRUE, ...)
        rtm <- rtm + bat.ret$rpt[1, 2] + td
        par <- bat.ret$par
        se2 <- bat.ret$se2
        phi <- par[1]
        ## mse <- knl.mnq.evl(rsp.dvp, knl.dvp, par, ...)[1, 2]
        mse <- knl.prd(rsp.dvp, knl.dvp, par, logged=FALSE, ...)[1, 2]
        ## record each iteration
        hst.num[[i+1]] <- c(ep=ep, bt=bt, rtm=rtm)
        hst.par[[i+1]] <- par
        hst.se2[[i+1]] <- se2
        hst.mse[[i+1]] <- mse

        ## message tracks
        if(i %% 80 == 0)
            cat(hd.msg(tks), "\n", sep="")
        cat(ln.msg(tks), "\n", sep="")

        if(i > max.itr) {cat('BMQ: reaching max iter:', max.itr, '\n'); break}
        if(rtm > wtm)   {cat('BMQ: reaching walltime:', wtm, 'h\n');    break}
        i <- i + 1
    }

    ## rearrange history
    hst.num <- do.call(rbind, hst.num)
    hst.mse <- do.call(rbind, hst.mse)
    hst.se2 <- DF(do.call(rbind, hst.se2))
    hst <- DF(hst.num, mse=hst.mse)
    
    ## mean parameter solution
    par <- Reduce(`+`, hst.par) / length(hst.par)

    ## mean SE
    hst.se2 <- split(hst.se2, hst.num[, 'ep'])
    hst.se2 <- sapply(hst.se2, function(.) 1 / colSums(1/.))
    se2 <- unname(rowMeans(hst.se2))

    rpt <- knl.prd(rsp.dvp, knl.dvp, par, logged=FALSE, ...)
    rpt <- rbind(DF(key='rtm', val=tail(hst$rtm, 1)), rpt)

    ## return the history and new parameters
    ret <- list(rpt=rpt, par=par, se2=se2, hst=hst)
    ret
}
