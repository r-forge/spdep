# Copyright 2002 by Roger Bivand 
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

LR.sarlm <- function(x, y)
{
	if (class(x) != "logLik") LLx <- logLik(x)
	else LLx <- c(x)
	if (class(y) != "logLik") LLy <- logLik(y)
	else LLy <- c(y)
	statistic <- 2*(LLx - LLy)
	attr(statistic, "names") <- "Likelihood ratio"
	parameter <- 1
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	estimate <- c(LLx, LLy)
	attr(estimate, "names") <- c(paste("Log likelihood of",
		deparse(substitute(x))), paste("Log likelihood of",
		deparse(substitute(y))))
	method <- "Likelihood ratio for spatial linear models"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res
}

logLik.sarlm <- function(object, ...) {
	LL <- c(object$LL)
	class(LL) <- "logLik"
	N <- length(residuals(object))
	attr(LL, "nall") <- N
	attr(LL, "nobs") <- N
	attr(LL, "df") <- object$parameters
	LL
}

