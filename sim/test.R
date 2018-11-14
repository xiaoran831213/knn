test.euc2 <- function(N=500, P=1000, t=20)
{
    library(microbenchmark)

    a <- matrix(rnorm(N * P), N, P)
    r <- microbenchmark(
        d1 <- as.matrix(dist(a))^2,
        d2 <- sqrt(abs(euc2(a)))^2,
        times=t)
    print(r)
    
    sum(abs(round(d1 - d2, 6)))
}


test.solve <- function(N=500, P=1000, t=20)
{
    library(microbenchmark)
    a <- tcrossprod(matrix(rnorm(N * P), N, P))
    b <- rnorm(N)

    inv <- function() {solve(a) %*% b}
    slv <- function() {solve(a, b)}
    qrs <- function() {qr.solve(a)}
    chl <- function() {chol2inv(chol(a)) %*% b}
    kpa <- function() {kappa(a)}

    r <- microbenchmark(inv, qrs, slv, chl, kpa, times=t)
    r
}

test.out3 <- function(N=500, P=1000, t=20)
{
    library(microbenchmark)
    a <- matrix(rnorm(N * P), N, P)

    f1 <- function()
    {
        (apply(a, 1L, function(x) sum(tcrossprod(x, x))) - rowSums(a^2))/2
    }
    
    f2 <- function()
    {
        r <- numeric(N)
        for(i in seq.int(P-1))
        {
            for(j in seq.int(i+1, P))
            {
                r <- r + a[, i] * a[, j]
            }
        }
        r
    }

    print(microbenchmark(r1 <- f1(), r2 <- f2(), times=t))
    
    list(r1, r2)
}

test.kpl1 <- function(N=500, P=2000)
{
    g <- readRDS('data/p35_c05.rds')$gmx
    x <- g[sample.int(nrow(g), N), sample.int(ncol(g), P)]
    x
}

source('sim/lnk.R')
plot.lnk <- function(d0, hist=0, seed=100)
{
    library(ggplot2)
    g <- ggplot()
    r <- 6
    g <- g + xlim(-r, +r) + ylim(-r, +r)
    g <- g + ylab(quote(bold(eta) %~% N(0, Sigma)))
    g <- g + xlab(quote(tilde(bold(eta))))

    ## histogram of normal distribution, plot in the back in gray
    if(hist == 1)
        g <- g + geom_histogram(aes(x=x, y=6*..density..), d0, fill='grey', alpha=1, binwidth=.5)

    g <- g + geom_point(aes(x=x, y=x), d0, color='grey', alpha=.01)
    g <- g + geom_point(aes(x=y, y=x, color=d), d0, alpha=.01)

    ## histogram of tranformed data, plot at the front in color
    if(hist == 1)
        g <- g + geom_histogram(aes(x=y, y=6*..density.., fill=d), d0, alpha=.3,
                                color='black', binwidth=.5)

    g <- g + coord_flip()
    g <- g + facet_wrap(~d, 1)
    g <- g + theme(
        legend.position = "none",
        strip.text=element_text(face='bold', size=15),
        strip.background = element_rect(colour="red", fill="#CCCCFF"),
        axis.title = element_text(face='bold', size=15))
    ## axis.title=element_blank()
    g
}

test.lnk <- function()
{
    x <- rnorm(10000, 0, sqrt(2.0))
    ## mapping by quantile
    dat <- within(list(),
    {
        ## st <- DF(d='st', x=x, y=dc(x, 'st', 0.05))
        ## ca <- DF(d='ca', x=x, y=dc(x, 'ca', 0.05))
        bn <- DF(d='binomial', x=x, y=dc(x, 'bn', 0.00))
        ch <- DF(d='chi-square', x=x, y=dc(x, 'ch', 0.00))
        ps <- DF(d='Poisson', x=x, y=dc(x, 'ps', 0.00))
        ex <- DF(d='exponential', x=x, y=dc(x, 'ex', 0.00))
        ## hy1 <- DF(d='Hyperbola 1',     x=x, y=HY1(x))
        ## hy2 <- DF(d='Hyperbola 2',     x=x, y=HY2(x))
        ## rc1 <- DF(d='Ricker Curver 1', x=x, y=RC1(x))
        ## rc2 <- DF(d='Ricker Curver 2', x=x, y=RC2(x))
    })
    d0 <- do.call(rbind, dat)
    ggsave('~/Dropbox/pub/img/mbq_tran.png', plot.lnk(d0, 0), scale=.85, width=19, height=19/4-.3)
    ggsave('~/Dropbox/pub/img/mbq_hist.png', plot.lnk(d0, 1), scale=.85, width=19, height=19/4-.3)

    ## non-linear
    dat <- within(list(),
    {
        hy1 <- DF(d='Hyperbola 1',     x=x, y=HY1(x))
        hy2 <- DF(d='Hyperbola 2',     x=x, y=HY2(x))
        rc1 <- DF(d='Ricker Curver 1', x=x, y=RC1(x))
        rc2 <- DF(d='Ricker Curver 2', x=x, y=RC2(x))
    })
    d0 <- do.call(rbind, dat)
    ggsave('~/Dropbox/pub/img/nlr_tran.png', plot.lnk(d0, 0), scale=.85, width=19, height=19/4-.3)
    ggsave('~/Dropbox/pub/img/nlr_hist.png', plot.lnk(d0, 1), scale=.85, width=19, height=19/4-.3)
}
