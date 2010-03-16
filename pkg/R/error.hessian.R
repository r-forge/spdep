getVmate_eig <- function(coefs, env,
#y, X, wy, WX, n, eig, 
    s2, trs, tol.solve=1.0e-10, optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_errlm_eig, env=env,
#y=y, X=X, wy=wy, WX=WX, n=n, eig=eig, 
            method="BFGS", hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_errlm_eig, env)
#, y, X, wy, WX, n, eig)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asye(coefs, env, s2, mat, trs)
#X, WX, n,             
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

f_errlm_eig <- function(coefs, env)
#y, X, wy, WX, n, eig) 
{
    lambda <- coefs[1]
    beta <- coefs[-1]
    yl <- get("y", envir=env) - lambda * get("wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    res <- yl - (xl %*% beta)
    SSE <- sum(res^2)
    n <- get("n", envir=env)
    s2 <- SSE/n
    eig <- get("eig", envir=env)
    if (is.complex(eig)) 
        det <- sum(log(1 - lambda * Re(eig)))
    else det <- sum(log(1 - lambda * eig))
    ret <- (det - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret,
        " Jacobian:", det, " SSE:", SSE, "\n")
   ret
}

insert_asye <- function(coefs, env, s2, mat, trs) {
#X, WX, n, s2, mat, trs) {
    lambda <- coefs[1]
    p <- length(coefs)-1
    p2 <- p+2
    omat <- matrix(0, nrow=p2, ncol=p2)
    LX <- get("x", envir=env) - lambda * get("WX", envir=env)
#X - lambda * WX
    omat[3:p2, 3:p2] <- -crossprod(LX)*s2
    omat[2, 2] <- mat[1, 1]
    n <- get("n", envir=env)
    omat[1, 1] <- -n/(2*(s2^2))
    omat[1, 2] <- omat[2, 1] <- -trB(lambda, trs)/s2
    omat
}

f_errlm_spam <- function(coefs, env)
#y, X, wy, WX, n, W, I) 
{
    lambda <- coefs[1]
    beta <- coefs[-1]
    yl <- get("y", envir=env) - lambda * get("wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    res <- yl - (xl %*% beta)
#    res <- (y - lambda * wy) - ((X - lambda * WX) %*% beta)
    SSE <- sum(res^2)
    n <- get("n", envir=env)
    s2 <- SSE/n
    W <- get("csrw", envir=env)
    I <- get("I", envir=env)
#    s2 <- SSE/n
    J1 <- try(determinant((I - lambda * W), logarithm=TRUE)$modulus,
        silent=TRUE)
   if (class(J1) == "try-error") {
      	Jacobian <- NA
   } else {
      	Jacobian <- J1
   }
   ret <- (Jacobian - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret,
        " Jacobian:", Jacobian, " SSE:", SSE, "\n")
   ret
}


getVmate_spam <- function(coefs, env,
#y, X, wy, WX, n, W, I,  
    s2, trs, tol.solve=1.0e-10, optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_errlm_spam, env=env,
#y=y, X=X, wy=wy, WX=WX, n=n, W=W, I=I, 
            method="BFGS", hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_errlm_spam, env=env)
#y, X, wy, WX, n, W, I)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asye(coefs, env, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

f_errlm_Matrix <- function(coefs, env)
#y, X, wy, WX, n, I, csrw, nW, nChol, pChol) 
{
    lambda <- coefs[1]
    beta <- coefs[-1]
    yl <- get("y", envir=env) - lambda * get("wy", envir=env)
    xl <- get("x", envir=env) - lambda * get("WX", envir=env)
    res <- yl - (xl %*% beta)
#    res <- (y - lambda * wy) - ((X - lambda * WX) %*% beta)
    SSE <- sum(res^2)
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
#    Jacobian <- determinant(I - lambda * csrw, logarithm=TRUE)$modulus
    ret <- (Jacobian - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret,
        " Jacobian:", Jacobian, " SSE:", SSE, "\n")
   ret
}

getVmate_Matrix <- function(coefs, env,
#y, X, wy, WX, n, I, csrw, nW, nChol, pChol,
    s2, trs, tol.solve=1.0e-10, optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_errlm_Matrix, env=env,
#y=y, X=X, wy=wy, WX=WX, n=n, I=I, csrw=csrw, nW=nW, nChol=nChol, pChol=pChol,
            method="BFGS", hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_errlm_Matrix, env)
#y, X, wy, WX, n, I, csrw, nW, nChol, pChol)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asye(coefs, env, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

