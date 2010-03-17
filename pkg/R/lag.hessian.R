f_laglm_eig <- function(coefs, env) {
    rho <- coefs[1]
    beta <- coefs[-1]
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml2_sse_env", env, rho, beta, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        res <- (get("y", envir=env) - rho * get("wy", envir=env)) - 
            get("x", envir=env) %*% beta
        SSE <- sum(res^2)
    }
    n <- get("n", envir=env)
    s2 <- SSE/n
    eig <- get("eig", envir=env)
    if (is.complex(eig)) 
        det <- sum(log(1 - rho * Re(eig)))
    else det <- sum(log(1 - rho * eig))
    ret <- (det - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("(eigen) rho:\t", rho, "\tfunction value:\t", ret, "\n")
   ret
}

getVmat_eig <- function(coefs, env,
    s2, trs, tol.solve=1.0e-10, optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_eig, env=env,
            method="BFGS", hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_laglm_eig, env)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asy(coefs, env,
            s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

f_laglm_Matrix <- function(coefs, env) {
    rho <- coefs[1]
    beta <- coefs[-1]
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml2_sse_env", env, rho, beta, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        res <- (get("y", envir=env) - rho * get("wy", envir=env)) - 
            get("x", envir=env) %*% beta
        SSE <- sum(res^2)
    }
    n <- get("n", envir=env)
    s2 <- SSE/n
    a <- -.Machine$double.eps^(1/2)
    b <- .Machine$double.eps^(1/2)
    .f <- if (package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1

    W <- get("W", envir=env)
    nW <- get("nW", envir=env)
    pChol <- get("pChol", envir=env)
    nChol <- get("nChol", envir=env)
    Jacobian <- ifelse(rho > b, n * log(rho) +
        (.f * c(determinant(update(nChol, nW, 1/rho))$modulus)),
        ifelse(rho < a, n* log(-(rho)) + (.f * c(determinant(update(pChol,
        W, 1/(-rho)))$modulus)), 0.0))
    ret <- (Jacobian - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("(Matrix) rho:\t", rho, "\tfunction value:\t", ret, "\n")
   ret
}

getVmat_Matrix <- function(coefs, env,
    s2, trs, tol.solve=1.0e-10, optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_Matrix, env=env,
            method="BFGS", hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_laglm_Matrix, env)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asy(coefs, env,
             s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

f_laglm_spam <- function(coefs, env) {
    rho <- coefs[1]
    beta <- coefs[-1]
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml2_sse_env", env, rho, beta, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        res <- (get("y", envir=env) - rho * get("wy", envir=env)) - 
            get("x", envir=env) %*% beta
        SSE <- sum(res^2)
    }
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
   ret <- (Jacobian - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
    if (get("verbose", envir=env)) cat("(spam) rho:\t", rho, "\tfunction value:\t", ret, "\n")
   ret
}


getVmat_spam <- function(coefs, env,
#y, X, yW, n, W, I, 
    s2, trs, tol.solve=1.0e-10, optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_spam, env=env,
#y=y, X=X, yW=yW, n=n,W=W, I=I, 
            method="BFGS", hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_laglm_spam, env)
#y, X, yW, n, W, I)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asy(coefs, env,
#y, X, yW, n, 
             s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

trB <- function(rho, tr)  sum(sapply(0:(length(tr)-1),
    function(i) rho^i * tr[i+1]))

insert_asy <- function(coefs, env,
#y, X, yW, n, 
    s2, mat, trs) {
    p <- length(coefs)-1
    p2 <- p+2
    n <- get("n", envir=env)
    omat <- matrix(0, nrow=p2, ncol=p2)
    omat[3:p2, 3:p2] <- -crossprod(get("x", envir=env))/s2
    omat[2, 2] <- mat[1, 1]
    omat[2, 3:p2] <- omat[3:p2, 2] <- -c(crossprod(get("wy", envir=env),
        get("x", envir=env))/s2)
    omat[1, 1] <- -n/(2*(s2^2))
    omat[1, 2] <- omat[2, 1] <- -trB(coefs[1], trs)/s2
    omat
}

