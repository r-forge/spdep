# Copyright 1998-2009 by Roger Bivand and Andrew Bernat
#

lagsarlm <- function(formula, data = list(), listw, 
	na.action, type="lag", method="eigen", quiet=NULL, 
	zero.policy=NULL, interval=c(-1,0.999), tol.solve=1.0e-10, 
	tol.opt=.Machine$double.eps^0.5, 
        fdHess=NULL, optimHess=FALSE, trs=NULL) {
        timings <- list()
        .ptime_start <- proc.time()
        if (is.null(quiet)) quiet <- !get("verbose", env = .spdepOptions)
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", env = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, 
		method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
        if (is.null(fdHess)) fdHess <- method != "eigen"
        stopifnot(is.logical(fdHess))
	can.sim <- as.logical(NA)
	if (listw$style %in% c("W", "S")) 
		can.sim <- spdep:::can.be.simmed(listw)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	switch(type, lag = if (!quiet) cat("\nSpatial lag model\n"),
	    mixed = if (!quiet) cat("\nSpatial mixed autoregressive model\n"),
	    stop("\nUnknown model type\n"))
	if (!quiet) cat("Jacobian calculated using ")
	switch(method, 
		eigen = if (!quiet) cat("neighbourhood matrix eigenvalues\n"),
	        Matrix = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		    stop("Matrix method requires symmetric weights")
		    if (listw$style %in% c("B", "C") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		    stop("Matrix method requires symmetric weights")
                    if (listw$style == "U") stop("U style not permitted, use C")
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
		stop("...\nUnknown method\n"))
	y <- model.extract(mf, "response")
	x <- model.matrix(mt, mf)
	if (NROW(x) != length(listw$neighbours))
		stop("Input data and weights have different dimensions")
	n <- NROW(x)
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	if (any(is.na(wy))) stop("NAs in lagged dependent variable")
	if (type != "lag") {
		# check if there are enough regressors
	        if (m > 1) {
			WX <- matrix(nrow=n,ncol=(m-(K-1)))
			for (k in K:m) {
				wx <- lag.listw(listw, x[,k], 
				    zero.policy=zero.policy)
				if (any(is.na(wx))) 
				    stop("NAs in lagged independent variable")
				WX[,(k-(K-1))] <- wx
			}
		}
		if (K == 2) {
         	    # unnormalized weight matrices
                	if (!(listw$style == "W")) {
 	      			intercept <- as.double(rep(1, n))
       	        		wx <- lag.listw(listw, intercept, 
					zero.policy = zero.policy)
                    		if (m > 1) {
                        		WX <- cbind(wx, WX)
                    		} else {
			      		WX <- matrix(wx, nrow = n, ncol = 1)
                    		}
                	} 
            	}   
		m1 <- m + 1
		mm <- NCOL(x) + NCOL(WX)
            	xxcolnames <- character(mm)
		for (k in 1:m) xxcolnames[k] <- xcolnames[k]
		for (k in m1:mm) 
		    xxcolnames[k] <- paste("lag.", xcolnames[k-mm+m], sep="")
		x <- cbind(x, WX)
		colnames(x) <- xxcolnames
		m <- NCOL(x)
		rm(wx, WX)
	}
# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
	LL_null_lm <- logLik(lm(y ~ 1))
	m <- NCOL(x)
	similar <- FALSE
	lm.null <- lm(y ~ x - 1)
	lm.w <- lm.fit(x, wy)
	e.null <- lm.null$residuals
	e.w <- lm.w$residuals
	e.a <- t(e.null) %*% e.null
	e.b <- t(e.w) %*% e.null
	e.c <- t(e.w) %*% e.w
        env <- new.env(parent=globalenv())
        assign("y", y, envir=env)
        assign("wy", wy, envir=env)
        assign("x", x, envir=env)
        assign("n", n, envir=env)
        assign("m", m, envir=env)
        assign("K", K, envir=env)
        assign("e.a", e.a, envir=env)
        assign("e.b", e.b, envir=env)
        assign("e.c", e.c, envir=env)
        assign("verbose", !quiet, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	if (method == "eigen") {
		if (!quiet) cat("Computing eigenvalues ...\n")
		if (listw$style %in% c("W", "S") && can.sim) {
#			eig <- eigenw(similar.listw(listw))
                        eig <- eigen(similar.listw_Matrix(listw),
                            only.values=TRUE)$values
			similar <- TRUE
		} else eig <- eigenw(listw)
		if (!quiet) cat("\n")
#range inverted 031031, email from Salvati Nicola (and Rein Halbersma)
		if (is.complex(eig)) eig.range <- 1/range(Re(eig))
		else eig.range <- 1/range(eig)
                interval <- c(eig.range[1]+.Machine$double.eps,
                    eig.range[2]-.Machine$double.eps)
                assign("eig", eig, envir=env)
                timings[["eigen_set_up"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
		opt <- optimize(sar.lag.mixed.f, 
			lower=eig.range[1]+.Machine$double.eps, 
			upper=eig.range[2]-.Machine$double.eps, maximum=TRUE,
			tol=tol.opt, env=env)
#eig=eig, e.a=e.a, e.b=e.b, e.c=e.c, n=n, quiet=quiet)
		rho <- opt$maximum
		names(rho) <- "rho"
		LL <- opt$objective
		optres <- opt
                timings[["eigen_opt"]] <- proc.time() - .ptime_start
	} else {
	    if (method == "spam") {
        	if (listw$style %in% c("W", "S") & can.sim) {
	    		W <- listw2U_spam(similar.listw_spam(listw))
	    		similar <- TRUE
		} else W <- as.spam.listw(listw)
        	I <- diag.spam(1, n, n)
                assign("W", W, envir=env)
                assign("I", I, envir=env)
                timings[["spam_set_up"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
	        opt <- optimize(sar.lag.mix.f.sp,
		    interval=interval, maximum=TRUE, tol=tol.opt, env=env)
	        rho <- c(opt$maximum)
	        names(rho) <- "rho"
 	        LL <- c(opt$objective)
	        optres <- opt
                timings[["spam_opt"]] <- proc.time() - .ptime_start
        } else if (method == "Matrix") {
        	if (listw$style %in% c("W", "S") & can.sim) {
	    		W <- listw2U_Matrix(similar.listw_Matrix(listw))
	    		similar <- TRUE
		} else W <- as_dsTMatrix_listw(listw)
		W <- as(W, "CsparseMatrix")
        	I <- as_dsCMatrix_I(n)
		Imult <- 2
		if (listw$style == "B") {
                    Imult <- ceiling((2/3)*max(apply(W, 1, sum)))
		    interval <- c(-0.5, +0.25)
		} else interval <- c(-1.2, +1)
                nW <- - W
		pChol <- Cholesky(W, super=FALSE, Imult = Imult)
		nChol <- Cholesky(nW, super=FALSE, Imult = Imult)
                assign("W", W, envir=env)
                assign("nW", nW, envir=env)
                assign("pChol", pChol, envir=env)
                assign("nChol", nChol, envir=env)
                timings[["Matrix_set_up"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
	        opt <- optimize(sar.lag.mix.f.M,
		    interval=interval, maximum=TRUE, tol=tol.opt, env=env)
	        rho <- c(opt$maximum)
	        names(rho) <- "rho"
 	        LL <- c(opt$objective)
	        optres <- opt
                timings[["Matrix_opt"]] <- proc.time() - .ptime_start
	    }
	}
        .ptime_start <- proc.time()
	lm.lag <- lm((y - rho*wy) ~ x - 1)
	r <- residuals(lm.lag)
	fit <- y - r
	names(r) <- names(fit)
	coef.rho <- coefficients(lm.lag)
	names(coef.rho) <- colnames(x)
	SSE <- deviance(lm.lag)
	s2 <- SSE/n
        timings[["coefs"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	if (method != "eigen") {
                coefs <- c(rho, coef.rho)
                if (fdHess && method == "Matrix") {
                    fdHess <- getVmat_Matrix(coefs, env,
#y, x, wy, n, W, I, nW, nChol, pChol, 
                        s2, trs, tol.solve=tol.solve, optim=optimHess)
                    if (is.null(trs)) {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("rho", colnames(x))
 		        rest.se <- sqrt(diag(fdHess)[-1])
		        rho.se <- sqrt(fdHess[1,1])
                    } else {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("sigma2", "rho", colnames(x))
 		        rest.se <- sqrt(diag(fdHess)[-c(1,2)])
		        rho.se <- sqrt(fdHess[2,2])
                    }
		    LMtest <- NULL
		    varb <- FALSE
		    ase <- FALSE
                    timings[["Matrix_fdHess"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                } else if (fdHess && method == "spam") {
                    fdHess <- getVmat_spam(coefs, env,
#y, x, wy, n, W, I, 
                        s2, trs, tol.solve=1.0e-10, optim=optimHess)
                    if (is.null(trs)) {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("rho", colnames(x))
 		        rest.se <- sqrt(diag(fdHess)[-1])
		        rho.se <- sqrt(fdHess[1,1])
                    } else {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("sigma2", "rho", colnames(x))
 		        rest.se <- sqrt(diag(fdHess)[-c(1,2)])
		        rho.se <- sqrt(fdHess[2,2])
                    }
		    LMtest <- NULL
		    varb <- FALSE
		    ase <- FALSE
                    timings[["spam_fdHess"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
               } else {
		    rest.se <- NULL
		    rho.se <- NULL
		    LMtest <- NULL
		    ase <- FALSE
		    varb <- FALSE
               }
	} else {
                if (fdHess) {
                    coefs <- c(rho, coef.rho)
                    fdHess <- getVmat_eig(coefs, env,
#y, x, wy, n, eig, 
                       s2, trs, tol.solve=tol.solve, optim=optimHess)
                    if (is.null(trs)) {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("rho", colnames(x))
                    } else {
                        rownames(fdHess) <- colnames(fdHess) <- 
                            c("sigma2", "rho", colnames(x))
                    }
                    timings[["eigen_fdHess"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
		tr <- function(A) sum(diag(A))
# beware of complex eigenvalues!
		O <- (eig/(1-rho*eig))^2
		omega <- sum(O)
		if (is.complex(omega)) omega <- Re(omega)
		W <- listw2mat(listw)
		A <- solve(diag(n) - rho*W)
		AW <- A %*% W
		zero <- rbind(rep(0,length(coef.rho)))
		xtawxb <- s2*(t(x) %*% AW %*% x %*% coef.rho)
		V <- s2*(s2*tr(t(AW) %*% AW) +
			t(AW %*% x %*% coef.rho) %*%
			(AW %*% x %*% coef.rho)) + omega*s2^2
		inf1 <- rbind(n/2, s2*tr(AW), t(zero))
		inf2 <- rbind(s2*tr(AW), V, xtawxb)
		xtx <- s2*t(x) %*% x
		inf3 <- rbind(zero, t(xtawxb), xtx)
		inf <- cbind(inf1, inf2, inf3)
		varb <- (s2^2) * solve(inf, tol=tol.solve)
		rownames(varb) <- colnames(varb) <- 
			c("sigma", "rho", colnames(x))
		rest.se <- sqrt(diag(varb))[-c(1:2)]
		rho.se <- sqrt(varb[2,2])
		TW <- (W %*% W) + (t(W) %*% W)
		T22 <- sum(diag(TW))
		T21A <- sum(diag(TW %*% A))
		LMtest <- ((t(r) %*% W %*% r)/s2)^2
		LMtest <- LMtest/(T22 - ((T21A^2)*(rho.se^2)))
		ase <- TRUE
                timings[["eigen_se"]] <- proc.time() - .ptime_start
	}
	call <- match.call()
	ret <- structure(list(type=type, rho=rho, 
		coefficients=coef.rho, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), lm.model=lm.null,
		method=method, call=call, residuals=r, opt=optres,
		lm.target=lm.lag, fitted.values=fit,
		se.fit=NULL, formula=formula, similar=similar,
		ase=ase, rho.se=rho.se, LMtest=LMtest, 
		resvar=varb, zero.policy=zero.policy, aliased=aliased,
                listw_style=listw$style, interval=interval, fdHess=fdHess,
                optimHess=optimHess, insert=!is.null(trs),
                LLNullLlm=LL_null_lm,
                timings=do.call("rbind", timings)[, c(1, 3)]), class=c("sarlm"))
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

sar.lag.mixed.f <- function(rho, env)
#eig, e.a, e.b, e.c, n, quiet)
{
        e.a <- get("e.a", envir=env)
        e.b <- get("e.b", envir=env)
        e.c <- get("e.c", envir=env)
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
        n <- get("n", envir=env)
        eig <- get("eig", envir=env)
	s2 <- SSE/n
	if (is.complex(eig)) det <- Re(prod(1 - rho*eig)) 
	else det <- prod(1 - rho*eig)
	ret <- (log(det) - ((n/2)*log(2*pi)) - (n/2)*log(s2)
		- (1/(2*s2))*SSE)
	if (get("verbose", envir=env)) cat("(eigen) rho:\t", rho, "\tfunction value:\t", ret, "\n")
	ret
}



sar.lag.mix.f.sp <- function(rho, env) {
#W, I, e.a, e.b, e.c, n, quiet) {
        e.a <- get("e.a", envir=env)
        e.b <- get("e.b", envir=env)
        e.c <- get("e.c", envir=env)
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
        n <- get("n", envir=env)
	s2 <- SSE/n
        W <- get("W", envir=env)
        I <- get("I", envir=env)
	J1 <- try(determinant((I - rho * W), logarithm=TRUE)$modulus,
            silent=TRUE)
        if (class(J1) == "try-error") {
        	Jacobian <- NA
        } else {
        	Jacobian <- J1
        }
	ret <- (Jacobian
		- ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
	if (get("verbose", envir=env)) 
	    cat("(spam) rho:\t", rho, "\tfunction value:\t", ret, "\n")
	ret
}

sar.lag.mix.f.M <- function(rho, env) {
#W, I, e.a, e.b, e.c, n, nW, nChol, 
#		pChol, quiet) {
        e.a <- get("e.a", envir=env)
        e.b <- get("e.b", envir=env)
        e.c <- get("e.c", envir=env)
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
        n <- get("n", envir=env)
	s2 <- SSE/n
        .f <- if (package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1
        W <- get("W", envir=env)
        nW <- get("nW", envir=env)
        pChol <- get("pChol", envir=env)
        nChol <- get("nChol", envir=env)
        if (isTRUE(all.equal(rho, 0))) {
            Jacobian <- rho
        } else if (rho > 0) {
	    detTRY <- try(determinant(update(nChol, nW, 1/rho))$modulus,
                silent=TRUE)
            if (class(detTRY) == "try-error") {
                Jacobian <- NaN
            } else {
                Jacobian <- n * log(rho) + (.f * detTRY)
            }
	} else {
            detTRY <- try(determinant(update(pChol, W, 1/(-rho)))$modulus,
                silent=TRUE)
            if (class(detTRY) == "try-error") {
               Jacobian <- NaN
            } else {
               Jacobian <- n * log(-(rho)) + (.f * detTRY)
            }
	}	
#	Jacobian <- determinant(I - rho * W, logarithm=TRUE)$modulus
	ret <- (Jacobian
		- ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
	if (get("verbose", envir=env)) 
	    cat("(Matrix) rho:\t", rho, "\tfunction value:\t", ret, "\n")
	ret
}


