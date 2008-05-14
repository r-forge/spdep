# Copyright 1998-2008 by Roger Bivand and Andrew Bernat
#

lagsarlm <- function(formula, data = list(), listw, 
	na.action=na.fail, type="lag", method="eigen", quiet=TRUE, 
	zero.policy=FALSE, interval=c(-1,0.999), tol.solve=1.0e-10, 
	tol.opt=.Machine$double.eps^0.5#, cholAlloc=NULL
	) {
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, 
		method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
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
	m <- NCOL(x)
	similar <- FALSE
	if (method == "eigen") {
		if (!quiet) cat("Computing eigenvalues ...\n")
		if (listw$style %in% c("W", "S") && can.sim) {
			eig <- eigenw(similar.listw(listw))
			similar <- TRUE
		} else eig <- eigenw(listw)
		if (!quiet) cat("\n")
#range inverted 031031, email from Salvati Nicola (and Rein Halbersma)
		if (is.complex(eig)) eig.range <- 1/range(Re(eig))
		else eig.range <- 1/range(eig)
		lm.null <- lm(y ~ x - 1)
		lm.w <- lm.fit(x, wy)
		e.null <- lm.null$residuals
		e.w <- lm.w$residuals
		e.a <- t(e.null) %*% e.null
		e.b <- t(e.w) %*% e.null
		e.c <- t(e.w) %*% e.w
		opt <- optimize(sar.lag.mixed.f, 
			lower=eig.range[1]+.Machine$double.eps, 
			upper=eig.range[2]-.Machine$double.eps, maximum=TRUE,
			tol=tol.opt, eig=eig, e.a=e.a, e.b=e.b, e.c=e.c,
			n=n, quiet=quiet)
		rho <- opt$maximum
		names(rho) <- "rho"
		LL <- opt$objective
		optres <- opt
	} else {
		opt <- dosparse(listw=listw, y=y, x=x, wy=wy, K=K, quiet=quiet,
			tol.opt=tol.opt, method=method, interval=interval, 
			can.sim=can.sim,  
			zero.policy=zero.policy)
		rho <- c(opt$maximum)
		names(rho) <- "rho"
		LL <- c(opt$objective)
		similar <- opt$similar
		optres <- opt$opt
	}
	lm.lag <- lm((y - rho*wy) ~ x - 1)
	r <- residuals(lm.lag)
	fit <- y - r
	names(r) <- names(fit)
	coef.rho <- coefficients(lm.lag)
	names(coef.rho) <- colnames(x)
	SSE <- deviance(lm.lag)
	s2 <- SSE/n
	if (method != "eigen") {
		LLs <- opt$LLs
		lm.null <- opt$lm.null
		rest.se <- NULL
		rho.se <- NULL
		LMtest <- NULL
		ase <- FALSE
		varb <- FALSE
	} else {
		LLs <- NULL
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
	}
	call <- match.call()
	ret <- structure(list(type=type, rho=rho, 
		coefficients=coef.rho, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), lm.model=lm.null,
		method=method, call=call, residuals=r, opt=optres,
		lm.target=lm.lag, fitted.values=fit,
		se.fit=NULL, formula=formula, similar=similar,
		ase=ase, LLs=LLs, rho.se=rho.se, LMtest=LMtest, 
		resvar=varb, zero.policy=zero.policy, aliased=aliased),
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

sar.lag.mixed.f <- function(rho, eig, e.a, e.b, e.c, n, quiet)
{
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
	s2 <- SSE/n
	if (is.complex(eig)) det <- Re(prod(1 - rho*eig)) 
	else det <- prod(1 - rho*eig)
	ret <- (log(det) - ((n/2)*log(2*pi)) - (n/2)*log(s2)
		- (1/(2*s2))*SSE)
	if (!quiet) cat("(eigen) rho:\t", rho, "\tfunction value:\t", ret, "\n")
	ret
}



sar.lag.mix.f.sp <- function(rho, W, I, e.a, e.b, e.c, n, quiet) {
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
	s2 <- SSE/n
	J1 <- try(determinant((I - rho * W), logarithm=TRUE)$modulus,
            silent=TRUE)
        if (class(J1) == "try-error") {
        	Jacobian <- NA
        } else {
        	Jacobian <- J1
        }
	ret <- (Jacobian
		- ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
	if (!quiet) 
	    cat("(spam) rho:\t", rho, "\tfunction value:\t", ret, "\n")
	ret
}

sar.lag.mix.f.M <- function(rho, W, I, e.a, e.b, e.c, n, nW, nChol, 
		pChol, quiet) {
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
	s2 <- SSE/n
        if (isTRUE(all.equal(rho, 0))) {
            Jacobian <- rho
        } else if (rho > 0) {
	    detTRY <- try(Matrix:::ldetL2up(nChol, nW, 1/rho),
                silent=TRUE)
            if (class(detTRY) == "try-error") {
                Jacobian <- NaN
            } else {
                Jacobian <- n * log(rho) + detTRY
            }
	} else {
            detTRY <- try(Matrix:::ldetL2up(pChol, W, 1/(-rho)),
                silent=TRUE)
            if (class(detTRY) == "try-error") {
               Jacobian <- NaN
            } else {
               Jacobian <- n * log(-(rho)) + detTRY
            }
	}	
#	Jacobian <- determinant(I - rho * W, logarithm=TRUE)$modulus
	ret <- (Jacobian
		- ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
	if (!quiet) 
	    cat("(Matrix) rho:\t", rho, "\tfunction value:\t", ret, "\n")
	ret
}

dosparse <- function (listw, y, x, wy, K, quiet, tol.opt, method, interval, 
	can.sim, zero.policy=FALSE) {
	similar <- FALSE
	m <- ncol(x)
	n <- nrow(x)
	if (method == "spam") {
        	if (listw$style %in% c("W", "S") & can.sim) {
	    		W <- as.spam.listw(listw2U(similar.listw(listw)))
	    		similar <- TRUE
		} else W <- as.spam.listw(listw)
        	I <- diag.spam(1, n, n)
	} else if (method == "Matrix") {
        	if (listw$style %in% c("W", "S") & can.sim) {
	    		W <- as_dsTMatrix_listw(listw2U(similar.listw(listw)))
	    		similar <- TRUE
		} else W <- as_dsTMatrix_listw(listw)
		W <- as(W, "CsparseMatrix")
        	I <- as_dsCMatrix_I(n)
		Imult <- 2
		if (listw$style == "B") {
                    Imult <- ceiling((2/3)*max(apply(W, 1, sum)))
		    interval <- c(-0.5, +0.25)
		} else interval <- c(-2, +1)
                nW <- - W
		pChol <- Cholesky(W, super=FALSE, Imult = Imult)
		nChol <- Cholesky(nW, super=FALSE, Imult = Imult)
		ns1 <- last <- 10
		prho1 <- seq(sqrt(.Machine$double.eps), interval[2],
                    length.out=ns1)
		
		while (last >= ns1) {
                   pdet1 <- Matrix:::ldetL2up(nChol, nW, 1/prho1)
		   wp1 <- which(is.finite(pdet1))
		   last <- wp1[length(wp1)]
		   if (last == ns1) prho1 <- seq(interval[2], 
		       1.5*interval[2], length.out=ns1)
		}
                lwp1n <- prho1[last]
                lwp2n <- prho1[last+1]
		prho2 <- seq(lwp2n, lwp1n, length.out=ns1)
                pdet2 <- Matrix:::ldetL2up(nChol, nW, 1/prho2)
		wp2 <- which(is.finite(pdet2))
                lwp2n <- prho2[wp2[length(wp2)]]
		
		nrho1 <- seq(interval[1], -sqrt(.Machine$double.eps),
                    length.out=ns1)
		
		first <- 1
		while (first == 1) {
                   ndet1 <- Matrix:::ldetL2up(pChol, W, 1/(-nrho1))
		   wn1 <- which(is.finite(ndet1))
		   first <- wn1[1]
		   if (first == 1) prho1 <- seq(1.5*interval[1], 
			interval[1], length.out=ns1)
		}

                lwn1n <- nrho1[wn1[1]]
                lwn2n <- nrho1[wn1[1]-1]
		nrho2 <- seq(lwn2n, lwn1n, length.out=ns1)
                ndet2 <- Matrix:::ldetL2up(pChol, W, 1/(-nrho2))
		wn2 <- which(is.finite(ndet2))
                lwn2n <- nrho2[wn2[1]]
		interval <- c(lwn2n, lwp2n)
		if (!quiet) cat("using interval:", interval, "\n")
	}
	LLs <- NULL
	# intercept-only bug fix Larry Layne 20060404
	if (m > 1) {
	    LLs <- vector(mode="list", length=length(K:m))
	    j <- 1
	    for (i in K:m) {
		# drop bug found by Gilles Spielvogel 20050128
		thisx <- x[,-i, drop = FALSE]
		lm.null <- lm.fit(thisx, y)
		lm.w <- lm.fit(thisx, wy)
		e.null <- lm.null$residuals
		e.w <- lm.w$residuals
		e.a <- t(e.null) %*% e.null
		e.b <- t(e.w) %*% e.null
		e.c <- t(e.w) %*% e.w
		if (method == "spam") {
		    LLs[[j]] <- optimize(sar.lag.mix.f.sp,
			interval=interval, maximum=TRUE, tol=tol.opt, W=W, I=I,
			e.a=e.a, e.b=e.b, e.c=e.c, n=n, quiet=quiet)$objective
		} else if (method == "Matrix") {
		    LLs[[j]] <- optimize(sar.lag.mix.f.M,
			interval=interval, maximum=TRUE, tol=tol.opt, W=W, I=I,
			e.a=e.a, e.b=e.b, e.c=e.c, n=n, nW=nW, nChol=nChol, 
			pChol=pChol, quiet=quiet)$objective
		}
#		gc(FALSE)
		attr(LLs[[j]], "nall") <- n
		attr(LLs[[j]], "nobs") <- n
		attr(LLs[[j]], "df") <- (m+2)-1
		attr(LLs[[j]], "name") <- colnames(x)[i]
		class(LLs[[j]]) <- "logLik"
		j <- j + 1
	    }
	}
	lm.null <- lm(y ~ x - 1)
	lm.w <- lm.fit(x, wy)
	e.null <- lm.null$residuals
	e.w <- lm.w$residuals
	e.a <- t(e.null) %*% e.null
	e.b <- t(e.w) %*% e.null
	e.c <- t(e.w) %*% e.w
	if (method == "spam") {
	    opt <- optimize(sar.lag.mix.f.sp,
		interval=interval, maximum=TRUE, tol=tol.opt, W=W, I=I,
		e.a=e.a, e.b=e.b, e.c=e.c, n=n, quiet=quiet)
	} else if (method == "Matrix") {
	    opt <- optimize(sar.lag.mix.f.M,
		interval=interval, maximum=TRUE, tol=tol.opt, W=W, I=I,
		e.a=e.a, e.b=e.b, e.c=e.c, n=n, nW=nW, nChol=nChol, 
		pChol=pChol, quiet=quiet)
	}
	maximum <- opt$maximum
	objective <- opt$objective
#	gc(FALSE)
	res <- list(maximum=maximum, objective=objective, LLs=LLs,
		lm.null=lm.null, similar=similar, opt=opt)
	res
}
