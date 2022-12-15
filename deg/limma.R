library(limma)
library(data.table)

fit <- function(d, s, f) {
    v <- voom(d, design=model.matrix(as.formula(f), data=s))
    l <- lmFit(v)
    l <- eBayes(l)
    r <- list(voom=v$E, wts=v$weights, coef=l$coef, t=l$t, p=l$p.value)
    dimnames(r$wts)<-dimnames(r$voom)
    l <- dimnames(r$coef)
    names(l) <- c(names(dimnames(d))[1], 'var')
    for(n in c('coef', 't', 'p'))
        dimnames(r[[n]]) <- l
    r
}

fit1 <- function(d, s, f) {
    l <- lmFit(d, design=model.matrix(as.formula(f), data=s))
    l <- eBayes(l)
    r <- list(coef=l$coef, t=l$t, p=l$p.value)
    l <- dimnames(r$coef)
    names(l) <- c(names(dimnames(d))[1], 'var')
    for(n in c('coef', 't', 'p'))
        dimnames(r[[n]]) <- l
    r
}

fit3 <- function(d, s, f) {
    s <- model.matrix(as.formula(f), data=s)
    v <- voom(d, design=s)
    l <- lmFit(v$E, design=s, weights=v$weights)
    r <- list(
        coef=l$coef,
        cov=l$cov,
        std=l$stdev.unscaled,
        sigma=matrix(l$sigma),
        df=matrix(l$df.residual)
    ) #, pivot=l$pivot)
    names(dimnames(r$coef)) <- c(names(dimnames(d))[1], 'var')
    names(dimnames(r$std)) <- names(dimnames(r$coef))
    names(dimnames(r$cov)) <- c('var', 'var1')
    dimnames(r$sigma) <- c(dimnames(d)[1], list('d0'=c('0')))
    dimnames(r$df) <- dimnames(r$sigma)
    r
}

fit4 <- function(coef, cov, std, sigma, df, contr) {
    l <- list(
        coefficients=coef,
        cov.coefficients=cov,
        stdev.unscaled=std,
        sigma=c(sigma),
        df.residual=c(df)
    )
    l <- contrasts.fit(l, contrasts=contr)
    l <- eBayes(l)
    r <- list(coef=l$coef, t=l$t, p=l$p.value)
    l <- dimnames(r$coef)
    names(l) <- c(names(dimnames(coef))[1], names(dimnames(contr))[2])
    for(n in c('coef', 't', 'p'))
        dimnames(r[[n]]) <- l
    r
}
