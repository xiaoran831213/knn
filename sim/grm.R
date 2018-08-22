#' read GRM binary of GCTA
readGRM <- function(pfx)
{
    ## get file names
    fn.rmx <- paste0(pfx, ".grm.bin")
    fn.N <- paste0(pfx, ".grm.N.bin")
    fn.id <- paste0(pfx, ".grm.id")

    ## 1) read families
    id <- read.table(fn.id, col.names=c('FID', 'IID'), as.is=TRUE)
    M <- nrow(id)                       # sample size
    L <- M * (M + 1) / 2                # number of GRs
    nms <- with(id, paste(FID, IID, sep='.'))
    
    ## unit size of a GR: 4 - single precision, 8 - double precision
    S <- file.size(fn.rmx) / L

    ## 2) read relatedness matrix
    rmx <- matrix(.0, M, M, dimnames=list(nms, nms))
    rmx[upper.tri(rmx, 1)] <- readBin(fn.rmx, .0, L, S)
    rmx[lower.tri(rmx, 0)] <- t(rmx)[lower.tri(rmx, 0)]

    ## 3) read variant count matrix
    S <- file.size(fn.N) / L
    if(S < 1)                           # same count for all
    {
        S <- file.size(fn.N)
        L <- 1
    }
    
    N <- matrix(0L, nrow(id), nrow(id), dimnames=list(nms, nms))
    N[upper.tri(N, 1)] <- readBin(fn.N, .0, L, S)
    N[lower.tri(N, 0)] <- t(N)[lower.tri(N, 0)]
    
    list(rmx=rmx, N=N, id=id)
}

saveGRM <- function(pfx, grm)
{
    ## get file names
    fn.rmx <- paste0(pfx, ".grm.bin")
    fn.N <- paste0(pfx, ".grm.N.bin")
    fn.id <- paste0(pfx, ".grm.id")

    ## complete id and N
    if(is.matrix(grm))
    {
        grm <- list(rmx=grm, id=makeID(grm), N=1.0)
    }

    with(grm,
    {
        ## upper.tri of col major = lower.tri of row major
        idx <- upper.tri(diag(nrow(id)), T)
        
        ## genomic relatedness matrix
        rmx <- rmx[idx]
        writeBin(rmx, fn.rmx, 4L)

        ## genomic variant count matrix
        N <- N[idx]
        writeBin(N, fn.N, 4L)

        ## subject IDs
        write(t(id), fn.id, 2, sep='\t')
    })
}

## mean of relationship
meanGRM <- function(x)
{
    ## the first GRM
    rmx <- .0
    N <- 0L
    id <- NA

    for(e in x)
    {
        if(is.character(e))
            e <- readGRM(e)
        rmx <- rmx + e$rmx * e$N
        N <- N + e$N
        id <- e$id
    }
    rmx <- rmx/N
    
    list(rmx=rmx, N=N, id=id)
}

makeID <- function(x)
{
    N <- nrow(x)
    fid <- sprintf('F%04X', seq(N))
    iid <- sprintf('I%04X', seq(N))
    data.frame(FID=fid, IID=iid, stringsAsFactors=FALSE)
}
