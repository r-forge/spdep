# Copyright 1998-2010 by Roger Bivand (non-W styles Rein Halbersma)
#

errorsarlm <- function(formula, data = list(), listw, na.action, 
	method="eigen", quiet=NULL, zero.policy=NULL, interval=c(-1,0.999), 
	tol.solve=1.0e-10, tol.opt=.Machine$double.eps^0.5,
        returnHcov=TRUE, pWOrder=250, fdHess=NULL,
        optimHess=FALSE, trs=NULL, LAPACK=FALSE, compiled_sse=FALSE) {
        timings <- list()
        .ptime_start <- proc.time()
        if (is.null(quiet)) quiet <- !get("verbose", env = .spdepOptions)
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", env = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
        if (is.null(fdHess)) fdHess <- method != "eigen"
        stopifnot(is.logical(fdHess))
        stopifnot(is.logical(LAPACK))
	can.sim <- as.logical(NA)
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	if (!quiet) cat(paste("\nSpatial autoregressive error model\n", 
		"Jacobian calculated using "))
	switch(method,
		eigen = if (!quiet) cat("neighbourhood matrix eigenvalues\n"),
	        Matrix = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		    stop("Matrix method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		    stop("Matrix method requires symmetric weights")
		    if (!quiet) cat("sparse matrix techniques using Matrix\n")
		},
	        spam = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		    stop("spam method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		    stop("spam method requires symmetric weights")
		    if (!quiet) cat("sparse matrix techniques using spam\n")
		},
		stop("...\n\nUnknown method\n"))
	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	n <- NROW(x)
	m <- NCOL(x)
# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
        LL_null_lm <- NULL
	if ("(Intercept)" %in% colnames(x)) LL_null_lm <- logLik(lm(y ~ 1))
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
# added no intercept Guillaume Blanchet 091103
	if (m > 1 || (m == 1 && K == 1)) {
	    WX <- matrix(nrow=n,ncol=(m-(K-1)))
	    for (k in K:m) {
		wx <- lag.listw(listw, x[,k], zero.policy=zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
	    }
	}
	if (K == 2) {
# modified to meet other styles, email from Rein Halbersma
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
		if (m > 1) WX <- cbind(wx, WX)
		else WX <- matrix(wx, nrow=n, ncol=1)
	}
	colnames(WX) <- xcolnames
	rm(wx)
	similar <- FALSE

        env <- new.env(parent=globalenv())
        assign("y", y, envir=env)
        assign("x", x, envir=env)
        assign("wy", wy, envir=env)
        assign("WX", WX, envir=env)
        assign("n", n, envir=env)
        assign("p", m, envir=env)
        assign("verbose", !quiet, envir=env)
        assign("compiled_sse", compiled_sse, envir=env)
        assign("first_time", TRUE, envir=env)
        assign("LAPACK", LAPACK, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	if (method == "eigen") {
		if (!quiet) cat("Computing eigenvalues ...\n")
		if (listw$style %in% c("W", "S") & can.sim) {
                        eig <- eigen(similar.listw_Matrix(listw),
                            only.values=TRUE)$values
			similar <- TRUE
		} else eig <- eigenw(listw)
		if (!quiet) cat("\n")
		if (is.complex(eig)) eig.range <- 1/range(Re(eig))
		else eig.range <- 1/range(eig)
                assign("eig", eig, envir=env)
                timings[["eigen_set_up"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                if (compiled_sse) {
                    ptr <- .Call("opt_error_init", PACKAGE="spdep")
#                    cfn <- function(ptr) .Call("opt_error_free",
#                        ptr, PACKAGE="spdep")
#                    reg.finalizer(ptr, cfn, onexit=FALSE)
                    assign("ptr", ptr, envir=env)
                }
		opt <- optimize(sar.error.f, 
			lower=eig.range[1]+.Machine$double.eps, 
			upper=eig.range[2]-.Machine$double.eps, maximum=TRUE,
			tol=tol.opt, env=env)
		lambda <- opt$maximum
		names(lambda) <- "lambda"
		LL <- opt$objective
                if (compiled_sse) {
                    .Call("opt_error_free", get("ptr", envir=env),
                        PACKAGE="spdep")
                }
                timings[["eigen_opt"]] <- proc.time() - .ptime_start
	} else if (method == "spam") {
        	if (listw$style %in% c("W", "S") & can.sim) {
	    		csrw <- listw2U_spam(similar.listw_spam(listw))
			    similar <- TRUE
		} else csrw <- as.spam.listw(listw)
                W <- as.spam.listw(listw)
        	I <- diag.spam(1, n, n)
                assign("csrw", csrw, envir=env)
                assign("I", I, envir=env)
                timings[["spam_set_up"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                if (compiled_sse) {
                    ptr <- .Call("opt_error_init", PACKAGE="spdep")
#                    cfn <- function(ptr) .Call("opt_error_free",
#                        ptr, PACKAGE="spdep")
#                    reg.finalizer(ptr, cfn, onexit=FALSE)
                    assign("ptr", ptr, envir=env)
                }
		opt <- optimize(sar.error.f.sp, interval=interval, 
			maximum=TRUE, tol=tol.opt, env=env)
		lambda <- opt$maximum
		names(lambda) <- "lambda"
		LL <- opt$objective
                if (compiled_sse) {
                    .Call("opt_error_free", get("ptr", envir=env),
                        PACKAGE="spdep")
                }
                timings[["spam_opt"]] <- proc.time() - .ptime_start
	} else if (method == "Matrix") {
        	if (listw$style %in% c("W", "S") & can.sim) {
	    	    csrw <- listw2U_Matrix(similar.listw_Matrix(listw))
	    	    similar <- TRUE
		} else csrw <- as_dsTMatrix_listw(listw)
		csrw <- as(csrw, "CsparseMatrix")
		Imult <- 2
		if (listw$style == "B") {
                    Imult <- ceiling((2/3)*max(apply(W, 1, sum)))
		    interval <- c(-0.5, +0.25)
		} else interval <- c(-1.2, +1)
                nW <- - csrw
		pChol <- Cholesky(csrw, super=FALSE, Imult = Imult)
		nChol <- Cholesky(nW, super=FALSE, Imult = Imult)
                assign("csrw", csrw, envir=env)
                assign("nW", nW, envir=env)
                assign("pChol", pChol, envir=env)
                assign("nChol", nChol, envir=env)
                W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        	I <- as_dsCMatrix_I(n)
                timings[["Matrix_set_up"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                if (compiled_sse) {
                    ptr <- .Call("opt_error_init", PACKAGE="spdep")
#                    cfn <- function(ptr) .Call("opt_error_free",
#                        ptr, PACKAGE="spdep")
#                    reg.finalizer(ptr, cfn, onexit=FALSE)
                    assign("ptr", ptr, envir=env)
                }
		opt <- optimize(sar.error.f.M, interval=interval, 
			maximum=TRUE, tol=tol.opt, env=env)
		lambda <- opt$maximum
		names(lambda) <- "lambda"
		LL <- opt$objective
                if (compiled_sse) {
                    .Call("opt_error_free", get("ptr", envir=env),
                        PACKAGE="spdep")
                }
                timings[["Matrix_opt"]] <- proc.time() - .ptime_start
	}
        .ptime_start <- proc.time()
	lm.target <- lm(I(y - lambda*wy) ~ I(x - lambda*WX) - 1)
	r <- as.vector(residuals(lm.target))
	fit <- as.vector(y - r)
	p <- lm.target$rank
	SSE <- deviance(lm.target)
	s2 <- SSE/n
	rest.se <- (summary(lm.target)$coefficients[,2])*sqrt((n-p)/n)
	coef.lambda <- coefficients(lm.target)
	names(coef.lambda) <- xcolnames
	lm.model <- lm(formula, data)
	ase <- FALSE
	lambda.se <- NULL
	LMtest <- NULL
	asyvar1 <- FALSE
        Hcov <- NULL
        timings[["coefs"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        assign("first_time", TRUE, envir=env)
	if (method == "eigen") {
		tr <- function(A) sum(diag(A))
		W <- listw2mat(listw)
		A <- solve(diag(n) - lambda*W)
		WA <- W %*% A
		asyvar <- matrix(0, nrow=2+p, ncol=2+p)
		asyvar[1,1] <- n / (2*(s2^2))
		asyvar[2,1] <- asyvar[1,2] <- tr(WA) / s2
		asyvar[2,2] <- tr(WA %*% WA) + tr(crossprod(WA))
# bug found 100224 German Muchnik Izon
#		asyvar[3:(p+2),3:(p+2)] <- s2*(t(x - lambda*WX) %*% 
                xl <- (x - lambda*WX)
		asyvar[3:(p+2),3:(p+2)] <- crossprod(xl)
		asyvar1 <- solve(asyvar, tol=tol.solve)
		rownames(asyvar1) <- colnames(asyvar1) <- 
			c("sigma", "lambda", xcolnames)
		
		lambda.se <- sqrt(asyvar1[2,2])
                timings[["eigen_se"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                if (returnHcov) {
                    pp <- lm.model$rank
                    p1 <- 1L:pp
                    R <- chol2inv(lm.model$qr$qr[p1, p1, drop = FALSE])
                    B <- tcrossprod(R, x) %*% A
                    A <- solve(diag(n) - lambda*t(W))
                    C <- A %*% x %*% R
                    Hcov <- B %*% C
                    attr(Hcov, "method") <- method
                    timings[["eigen_hcov"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
                if (fdHess) {
                    coefs <- c(lambda, coef.lambda)
                    if (compiled_sse) {
                        ptr <- .Call("hess_error_init", PACKAGE="spdep")
#                        cfn <- function(ptr) .Call("hess_error_free",
#                            ptr, PACKAGE="spdep")
#                        reg.finalizer(ptr, cfn, onexit=FALSE)
                        assign("ptr", ptr, envir=env)
                    }
                    fdHess <- getVmate_eig(coefs, env,
                       s2, trs, tol.solve=tol.solve, optim=optimHess)
                    if (compiled_sse) {
                        .Call("hess_error_free", get("ptr", envir=env),
                            PACKAGE="spdep")
                    }
                    if (is.null(trs)) {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("lambda", colnames(x))
                    } else {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("sigma2", "lambda", colnames(x))
                    }
                    timings[["eigen_fdHess"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
		ase <- TRUE
	} else {
                if (fdHess && method == "Matrix") {
                    coefs <- c(lambda, coef.lambda)
                    if (compiled_sse) {
                        ptr <- .Call("hess_error_init", PACKAGE="spdep")
#                        cfn <- function(ptr) .Call("hess_error_free",
#                            ptr, PACKAGE="spdep")
#                        reg.finalizer(ptr, cfn, onexit=FALSE)
                        assign("ptr", ptr, envir=env)
                    }
                    fdHess <- getVmate_Matrix(coefs, env,
                        s2, trs, tol.solve=tol.solve, optim=optimHess)
                    if (compiled_sse) {
                        .Call("hess_error_free", get("ptr", envir=env),
                            PACKAGE="spdep")
                    }
                    if (is.null(trs)) {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("lambda", colnames(x))
 		        rest.se <- sqrt(diag(fdHess)[-1])
		        lambda.se <- sqrt(fdHess[1,1])
                    } else {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("sigma2", "lambda", colnames(x))
 		        rest.se <- sqrt(diag(fdHess)[-c(1,2)])
		        lambda.se <- sqrt(fdHess[2,2])
                    }
                    timings[["Matrix_fdHess"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                } else if (fdHess && method == "spam") {
                    coefs <- c(lambda, coef.lambda)
                    if (compiled_sse) {
                        ptr <- .Call("hess_error_init", PACKAGE="spdep")
#                        cfn <- function(ptr) .Call("hess_error_free",
#                            ptr, PACKAGE="spdep")
#                        reg.finalizer(ptr, cfn, onexit=FALSE)
                        assign("ptr", ptr, envir=env)
                    }
                    fdHess <- getVmate_spam(coefs, env,
                        s2, trs, tol.solve=tol.solve, optim=optimHess)
                    if (compiled_sse) {
                        .Call("hess_error_free", get("ptr", envir=env),
                            PACKAGE="spdep")
                    }
                    if (is.null(trs)) {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("lambda", colnames(x))
 		        rest.se <- sqrt(diag(fdHess)[-1])
		        lambda.se <- sqrt(fdHess[1,1])
                    } else {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("sigma2", "lambda", colnames(x))
 		        rest.se <- sqrt(diag(fdHess)[-c(1,2)])
		        lambda.se <- sqrt(fdHess[2,2])
                    }
                    timings[["spam_fdHess"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }                    
                if (returnHcov) {
                    pp <- lm.model$rank
                    p1 <- 1L:pp
                    R <- chol2inv(lm.model$qr$qr[p1, p1, drop = FALSE])
                    B <- tcrossprod(R, x)
                    B1 <- as(powerWeights(W=W, rho=lambda, order=pWOrder,
                        X=B, tol=tol.solve), "matrix")
                    C <- x %*% R
                    C1 <- as(powerWeights(W=t(W), rho=lambda, order=pWOrder,
                        X=C, tol=tol.solve), "matrix")
                    Hcov <- B1 %*% C1
                    attr(Hcov, "method") <- method
                    timings[["sparse_hcov"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
        }
	call <- match.call()
	names(r) <- names(y)
	names(fit) <- names(y)
	ret <- structure(list(type="error", lambda=lambda,
		coefficients=coef.lambda, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), lm.model=lm.model, 
		method=method, call=call, residuals=r, lm.target=lm.target,
		opt=opt, fitted.values=fit, ase=ase, formula=formula,
		se.fit=NULL, resvar=asyvar1, similar=similar,
		lambda.se=lambda.se, LMtest=LMtest, zero.policy=zero.policy, 
		aliased=aliased, LLNullLlm=LL_null_lm, Hcov=Hcov,
                interval=interval, fdHess=fdHess,
                optimHess=optimHess, insert=!is.null(trs),
                timings=do.call("rbind", timings)[, c(1, 3)]),
                class=c("sarlm"))
	if (zero.policy) {
		zero.regs <- attr(listw$neighbours, 
			"region.id")[which(card(listw$neighbours) == 0)]
		if (length(zero.regs) > 0)
			attr(ret, "zero.regs") <- zero.regs
	}
	if (!is.null(na.act))
		ret$na.action <- na.act
	ret
}

sar.error.f <- function(lambda, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml_sse_env", env, lambda, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        yl <- get("y", envir=env) - lambda * get("wy", envir=env)
        xl <- get("x", envir=env) - lambda * get("WX", envir=env)
	xl.q <- qr.Q(qr(xl, LAPACK=get("LAPACK", envir=env)))
	xl.q.yl <- crossprod(xl.q, yl)
	SSE <- crossprod(yl) - crossprod(xl.q.yl)
    }
    n <- get("n", envir=env)
    s2 <- SSE/n
    eig <- get("eig", envir=env)
    if (is.complex(eig)) 
        det <- sum(log(1 - lambda * Re(eig)))
    else det <- sum(log(1 - lambda * eig))
    ret <- (det - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", det, " SSE:", SSE, "\n")
    ret
}

sar.error.f.sp <- function(lambda, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml_sse_env", env, lambda, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        yl <- get("y", envir=env) - lambda * get("wy", envir=env)
        xl <- get("x", envir=env) - lambda * get("WX", envir=env)
	xl.q <- qr.Q(qr(xl, LAPACK=get("LAPACK", envir=env)))
	xl.q.yl <- crossprod(xl.q, yl)
	SSE <- crossprod(yl) - crossprod(xl.q.yl)
    }
    n <- get("n", envir=env)
    s2 <- SSE/n
    csrw <- get("csrw", envir=env)
    I <- get("I", envir=env)
    J1 <- try(determinant((I - lambda * csrw), logarithm=TRUE)$modulus,
        silent=TRUE)
    if (class(J1) == "try-error") {
        Jacobian <- NA
    } else {
        Jacobian <- J1
    }
    ret <- (Jacobian -
	((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
	if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", Jacobian, " SSE:", SSE, "\n")
	ret
}

sar.error.f.M <- function(lambda, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml_sse_env", env, lambda, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        yl <- get("y", envir=env) - lambda * get("wy", envir=env)
        xl <- get("x", envir=env) - lambda * get("WX", envir=env)
	xl.q <- qr.Q(qr(xl, LAPACK=get("LAPACK", envir=env)))
	xl.q.yl <- crossprod(xl.q, yl)
	SSE <- crossprod(yl) - crossprod(xl.q.yl)
    }
    n <- get("n", envir=env)
    s2 <- SSE/n
    csrw <- get("csrw", envir=env)
    nW <- get("nW", envir=env)
    pChol <- get("pChol", envir=env)
    nChol <- get("nChol", envir=env)
    a <- -.Machine$double.eps^(1/2)
    b <- .Machine$double.eps^(1/2)
    .f <- if (package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1

    Jacobian <- ifelse(lambda > b, n * log(lambda) +
            (.f * c(determinant(update(nChol, nW, 1/lambda))$modulus)),
            ifelse(lambda < a, n* log(-(lambda)) + 
            (.f * c(determinant(update(pChol, csrw, 1/(-lambda)))$modulus)),
            0.0))

    ret <- (Jacobian -
		((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", Jacobian, " SSE:", SSE, "\n")
	ret
}

