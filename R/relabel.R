## as.mcmc.relabel <- function(x) as.mcmc(x$nuLam)

relabel <- function(obj, ...) UseMethod("relabel")

relabel.default <- function(obj, ...) stop("No default method.  Sorry.")

## relabel.cfabugs <- function(obj, random=FALSE, max.iter=100, ...){
##     n.chains <- obj$n.chains
##     sims <- obj$sims.array
##     draws <- NULL
##     for (i in seq(n.chains)) draws <- rbind(draws, sims[, i, ])
##     colnames(draws) <- dimnames(sims)[[3]]
##     rlbl <- relabel(draws, random=random, max.iter=max.iter)
##     return(rlbl)
## }

## relabel.bcfa <- function(obj, random=FALSE, max.iter=100, ...){
##     parnames <- grep("Lam", colnames(obj$draws[[1]]), value=TRUE)
##     obj$nu <- vector("list", length=length(obj$draws))
##     for (i in 1:length(obj$draws)){
##         r <- relabel(obj$draws[[i]][, parnames], random=random, max.iter=max.iter)
##         obj$draws[[i]][, parnames] <- r$nuLam
##         obj$nu[[i]] <- r$nu
##     }
##     obj$factor.idx <- r$factor.idx
##     return(obj)
## }

relabel.matrix <- function(obj, random=FALSE, max.iter=100, ...){
    ## Multivariate version
    ## require(mvtnorm)
    nms <- dimnames(obj)[[2]]
    n <- nrow(obj)
    factor.idx <- do.call("rbind", lapply(strsplit(nms, ","), function(x) as.numeric(gsub("[^[:digit:]]", "", x))))[,2]
    q <- max(factor.idx)
    col.idx <- lapply(1:q, function(j) which(factor.idx==j))
    if (random){
        nu <- matrix(2*rbinom(q*nrow(obj), 1, 0.5) - 1, ncol=q, nrow=nrow(obj))
    } else {
        nu <- do.call(rbind,
                      lapply(split(obj, row(obj)),
                             function(x, idx) sapply(idx, function(cols) sign(sum(x[cols]))),
                             idx=col.idx))
          ## nu <- t(apply(obj, 1, function(x, idx) sapply(idx, function(cols) sign(sum(x[cols]))), idx=col.idx))
    }
    nu.init <- nu
    m <- lapply(1:q, function(k, nu, obj, col.idx) colMeans(nu[, k]*obj[, col.idx[[k]]]),
                nu=nu, obj=obj, col.idx=col.idx)
    s <- lapply(1:q, function(k, nu, obj, col.idx) apply(nu[, k]*obj[, col.idx[[k]]], 2, sd),
                nu=nu, obj=obj, col.idx=col.idx)
    logdnorm <- logdnorm.p <- matrix(nrow=n, ncol=q)
    for (k in 1:q){
        cat("Initializing nu[", k, "]\n", sep="")
        ci <- col.idx[[k]]
        for (i in 1:n){
            logdnorm[i, k] <- -sum(dnorm(nu[i, k]*obj[i, ci], m[[k]], s[[k]], log=TRUE))
        }
    }
    not.converged <- TRUE
    iter <- 0

    ##==================================================
    ## Main iterative loop
    while (not.converged && (iter<max.iter)){
        iter <- iter + 1
        nu.old <- nu

        ## Update nu parameters
        for (k in 1:q){
            cat("Iteration ", iter, ": Updating nu[", k, "]\n", sep="")
            ci <- col.idx[[k]]
            for (i in 1:n){
                logdnorm.p[i, k] <- -sum(dnorm(-nu[i, k]*obj[i, ci], m[[k]], s[[k]], log=TRUE))
            }
            switch.sign <- logdnorm[, k] > logdnorm.p[, k]
            nu[switch.sign, k] <- -nu[switch.sign, k]
            logdnorm[switch.sign, k] <- logdnorm.p[switch.sign, k]
        }

        ## Update means and variances
        m <- lapply(1:q, function(k, nu, obj, col.idx) colMeans(nu[, k]*obj[, col.idx[[k]]]), nu=nu, obj=obj, col.idx=col.idx)
        s <- lapply(1:q, function(k, nu, obj, col.idx) apply(nu[, k]*obj[, col.idx[[k]]], 2, sd), nu=nu, obj=obj, col.idx=col.idx)

        if (sum(abs(nu.old - nu))==0) not.converged <- !not.converged
    }

    cat("\n")
    if (max.iter==iter) warning("Maximum number of iterations reached without convergence.")

    out <- list()
    out$nuLam <- nu[, factor.idx]*obj
    colnames(out$nuLam) <- colnames(obj)
    out$nu.init <- nu.init
    out$nu <- nu
    out$m <- m
    out$s <- s
    out$factor.idx <- factor.idx
    out$col.idx <- col.idx
    out$iter <- iter
    out$loss <- sum(logdnorm)
    out$converged <- !not.converged
    class(out) <- c("relabel")
    return(out)
}

