# Copyright 2006-7 by Roger Bivand
#

as_dgRMatrix_listw <- function(listw) {
	if(!inherits(listw, "listw")) stop("not a listw object")
	n <- length(listw$neighbours)
	cardw <- card(listw$neighbours)
	p0 <- as.integer(c(0, cumsum(cardw)))
	scard <- sum(cardw)
	z <- .Call("listw2dgR", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard), PACKAGE="spdep")
	res <- new("dgRMatrix", j=z[[1]], p=p0, Dim=as.integer(c(n, n)),
		x=z[[2]])
	res
}

as_dsTMatrix_listw <- function(listw) {
	if (!inherits(listw, "listw")) stop("not a listw object")
	if (!is.symmetric.glist(listw$neighbours, listw$weights))
		stop("not a symmetric matrix")
	n <- length(listw$neighbours)
	cardw <- card(listw$neighbours)
	scard <- sum(cardw)
	if (scard %% 2 != 0) stop("odd neighbours sum")
	z <- .Call("listw2dsT", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard/2), PACKAGE="spdep")

	res <- new("dsTMatrix", i=z[[1]], j=z[[2]], Dim=as.integer(c(n, n)),
		x=z[[3]])
	res
}

as_dsCMatrix_I <- function(n) {
	if (n < 1) stop("matrix must have positive dimensions")
	as(as(Diagonal(n), "symmetricMatrix"), "CsparseMatrix")
}

as_dsCMatrix_IrW <- function(W, rho) {
	stopifnot(is(W, "symmetricMatrix"))
	n <- dim(W)[1]
	as_dsCMatrix_I(n) - rho * W
}

Jacobian_W <- function(W, rho) {
	sum(2*log(diag(chol(as_dsCMatrix_IrW(W, rho)))))
}

# ww <- as(as_dgRMatrix_listw(w), "CsparseMatrix")
# zz <- 0.5 * (ww + t(ww))

# w1 <- nb2listw(eire.nb, style="B")
# d <- Diagonal(x=1/sqrt(card(eire.nb)))
# d1 <- as(as(d, "symmetricMatrix"), "CsparseMatrix")
# ww1 <- as(as_dgRMatrix_listw(w1), "CsparseMatrix")
# sim1 <- d1 %*% ww1 %*% d1
# range(eigen(sim1)$values)
# range(eigenw(w))

