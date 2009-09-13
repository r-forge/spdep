f_laglm_eig <- function(coefs, y, X, yW, n, eig) {
    rho <- coefs[1]
    beta <- coefs[-1]
    res <- (y - rho * yW) - X %*% beta
    SSE <- sum(res^2)
    s2 <- SSE/n
    if (is.complex(eig)) 
        det <- sum(log(1 - rho * Re(eig)))
    else det <- sum(log(1 - rho * eig))
    ret <- (det - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
   ret
}

getVmat_eig <- function(coefs, y, X, yW, n, eig, tol.solve=1.0e-10,
    optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_eig, y=y, X=X, yW=yW, n=n, eig=eig,
            method="BFGS", hessian=TRUE)
        res <- solve(-(opt$hessian), tol.solve=tol.solve)
    } else {
        fd <- fdHess(coefs, f_laglm_eig, y, X, yW, n, eig)
        res <- solve(-(fd$Hessian), tol.solve=tol.solve)
    }
    res
}

f_laglm_Matrix <- function(coefs, y, X, yW, n, W, I, nW, nChol, pChol) {
    rho <- coefs[1]
    beta <- coefs[-1]
    res <- (y - rho * yW) - X %*% beta
    SSE <- sum(res^2)
    s2 <- SSE/n
    a <- -.Machine$double.eps^(1/2)
    b <- .Machine$double.eps^(1/2)
    .f <- if (package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1

    Jacobian <- ifelse(rho > b, n * log(rho) +
        (.f * c(determinant(update(nChol, nW, 1/rho))$modulus)),
        ifelse(rho < a, n* log(-(rho)) + (.f * c(determinant(update(pChol,
        W, 1/(-rho)))$modulus)), 0.0))
    ret <- (Jacobian - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
   ret
}

getVmat_Matrix <- function(coefs, y, X, yW, n, W, I, nW, nChol, pChol,
    tol.solve=1.0e-10, optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_Matrix, y=y, X=X, yW=yW, n=n,
            W=W, I=I, nW=nW, nChol=nChol, pChol=pChol, method="BFGS",
            hessian=TRUE)
        res <- solve(-(opt$hessian), tol.solve=tol.solve)
    } else {
        fd <- fdHess(coefs, f_laglm_Matrix, y, X, yW, n, W, I, nW, nChol,
            pChol)
        res <- solve(-(fd$Hessian), tol.solve=tol.solve)
    }
    res
}

f_laglm_spam <- function(coefs, y, X, yW, n, W, I) {
    rho <- coefs[1]
    beta <- coefs[-1]
    res <- (y - rho * yW) - X %*% beta
    SSE <- sum(res^2)
    s2 <- SSE/n
    J1 <- try(determinant((I - rho * W), logarithm=TRUE)$modulus,
        silent=TRUE)
   if (class(J1) == "try-error") {
      	Jacobian <- NA
   } else {
      	Jacobian <- J1
   }
   ret <- (Jacobian - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
   ret
}


getVmat_spam <- function(coefs, y, X, yW, n, W, I, tol.solve=1.0e-10,
    optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_spam, y=y, X=X, yW=yW, n=n,
            W=W, I=I, method="BFGS", hessian=TRUE)
        res <- solve(-(opt$hessian), tol.solve=tol.solve)
    } else {
        fd <- fdHess(coefs, f_laglm_spam, y, X, yW, n, W, I)
        res <- solve(-(fd$Hessian), tol.solve=tol.solve)
    }
    res
}

