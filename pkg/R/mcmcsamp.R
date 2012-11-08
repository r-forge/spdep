MCMCsamp <- function(object, mcmc = 1L, verbose = NULL, ...) UseMethod("MCMCsamp")
# from lme4/R/AllGeneric.R
MCMCsamp.spautolm <- function(object, mcmc = 1L, verbose = NULL, ...,
    burnin=0L, thin=1L, seed=NA, verbose_step=100L, listw, control=list()) {
    con <- list(Imult=2, cheb_q=5, MC_p=16, MC_m=30, super=NULL,
        spamPivot="MMD", in_coef=0.1, type="MC",
        correct=TRUE, trunc=TRUE, SE_method="LU", nrho=200,
        interpn=2000, small_asy=TRUE, small=1500, SElndet=NULL,
        LU_order=FALSE)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))
    if (verbose) stopifnot(verbose_step > 0)
    if (!inherits(listw, "listw")) 
        stop("No neighbourhood list")
    method <- object$method
    family <- object$family
    if (family == "SMA" && method != "eigen") stop("SMA only for eigen method")
    X <- object$X
    N <- nrow(X)
    if (N != length(listw$neighbours))
	 stop("Input data and neighbourhood list have different dimensions")
    stopifnot(ncol(X) == length(object$fit$coefficients))
    weights <- object$weights
    stopifnot(length(weights) == N)
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) 
	can.sim <- can.be.simmed(listw)
    sum_lw <- sum(log(weights))
    env <- new.env()
    assign("Y", object$Y, envir=env)
    assign("X", X, envir=env)
    assign("n", N, envir=env)
    assign("weights", weights, envir=env)
    assign("can.sim", can.sim, envir=env)
    assign("family", family, envir=env)
    assign("method", method, envir=env)
    assign("verbose", verbose, envir=env)
    assign("listw", listw, envir=env)
    assign("sum_lw", sum_lw, envir=env)
    W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
    if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
        warning("Non-symmetric spatial weights in CAR model")
    assign("W", W, envir=env)
    I <- as_dsCMatrix_I(N)
    assign("I", I, envir=env)
    Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
        "CsparseMatrix")
    assign("Sweights", Sweights, envir=env)

    if (verbose) cat(paste("\nJacobian calculated using "))

    interval <- jacobianSetup(method, env, con, trs=object$trs,
        interval=object$interval)
    assign("interval", interval, envir=env)

    start <- c(object$lambda, object$fit$coefficients)
    V <- object$fdHess
    stopifnot(nrow(V) == ncol(V))
    stopifnot(nrow(V) == length(start))
    res <- MCMCmetrop1R(fun=f_spautolm_hess, theta.init=start, burnin=burnin,
        mcmc=mcmc, thin=thin, verbose=ifelse(verbose, verbose_step, 0),
        seed=seed, logfun=TRUE, V=V, env=env)
    colnames(res) <- c(names(object$lambda), names(object$fit$coefficients))
    res
}
