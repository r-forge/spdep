# Copyright 2001-2 by Roger Bivand 
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#

lm.LMtests <- function(model, listw, zero.policy=FALSE, test="LMerr") {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(class(model) != "lm") stop(paste(deparse(substitute(model)),
		"not an lm object"))
	N <- length(listw$neighbours)
	u <- as.vector(resid(model))
	if (N != length(u)) stop("objects of different length")
	if (is.null(attr(listw$weights, "W")) || !attr(listw$weights, "W"))
		warning("Spatial weights matrix not row standardized")

	tracew <- function (listw) {
		dlmtr <- 0
		n <- length(listw$neighbours)
		for (i in 1:n) {
			dij <- listw$neighbours[[i]]
			ndij <- length(dij)
			wdij <- listw$weights[[i]]
			for (j in 1:ndij) {
				k <- dij[j]
				if (k > i) {
				    dk <- which(listw$neighbours[[k]] == i)
				    if (dk > 0 &&
					dk <= length(listw$neighbours[[k]]))
					wdk <- listw$weights[[k]][dk]
					else wdk <- 0
					dlmtr <- dlmtr + (wdk * wdk) + 2 *
					(wdij[j] * wdk) + (wdij[j] * wdij[j])
				}
			}
		}
		dlmtr
	}

	y <- model.response(model.frame(model))
	X <- model.matrix(terms(model), model.frame(model))
	yhat <- as.vector(fitted(model))
	p <- model$rank
	p1 <- 1:p
	XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
	sigma2 <- (t(u) %*% u) / N
	T <- tracew(listw)
	Wu <- lag.listw(listw, u, zero.policy)
	Wy <- lag.listw(listw, y, zero.policy)
	Wyhat <- lag.listw(listw, yhat, zero.policy)
	XtWyhat <- t(X) %*% Wyhat
	dutWu <- (t(u) %*% Wu) / sigma2
	resa <- (dutWu ^ 2) / T
	J <- (1/(N*sigma2)) *
		((t(Wyhat) %*% Wyhat) -
		(t(XtWyhat) %*% XtXinv %*% XtWyhat) +
		(T * sigma2))
	dutWy <- (t(u) %*% Wy) / sigma2
	nt <- length(test)
	tres <- vector(mode="list", length=nt)
	names(tres) <- test
	for (i in 1:nt) {
		testi <- test[i]
		zz <- switch(testi,
		LMerr = vec <- c(resa, 1),
		LMlag = vec <- c((dutWy ^ 2) / (N * J), 1),
		RLMerr = vec <- c(((dutWu - (T*((N*J)^-1))*dutWy)^2) /
			(T * (1 - T*((N*J)^-1))), 1),
		RLMlag = vec <- c(((dutWy - dutWu)^2)/ ((N*J) - T), 1),
		SARMA = vec <- c(((dutWy - dutWu)^2)/ ((N*J) - T) + resa, 2)
		)
		if (is.null(zz)) stop(paste(testi, ": no such test", sep=""))
		statistic <- vec[1]
		names(statistic) <- testi
		parameter <- vec[2]
		names(parameter) <- "df"
		p.value <- 1 - pchisq(statistic, parameter)
		names(p.value) <- ""
		method <- "Lagrange multiplier diagnostics for spatial dependence"
		data.name <- paste("model:", deparse(model$call),
    	    	"\nweights:", deparse(substitute(listw)), "\n")
		tres[[i]] <- list(statistic=statistic, parameter=parameter,
			p.value=p.value, method=method, data.name=data.name)
		class(tres[[i]]) <- "htest"
	}
	class(tres) <- "LMtestlist"
	tres
}

print.LMtestlist <- function(x, ...) {
	for (i in 1:length(x)) print(x[[i]])
	invisible(x)
}
