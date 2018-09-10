source('R/gsm.R')

test.gsm <- function(N=9L, P=6L, ...)
{
    ## genotype
    G <- matrix(sample.int(3L, P * N, TRUE) - 1L, N, P)
    colnames(G) <- sprintf('%02d', 1L:P)

    ## model
    M <- ~ a*a[]

    ## expanded input matrix
    X <- gsm(M, G=G, rm.nic=TRUE, ...)

    ## return
    list(G=G, X=X, M=mdl.str(M))
}
