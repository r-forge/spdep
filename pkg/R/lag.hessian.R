f_laglm_eig <- function(coefs, n, eig) {
    rho <- coefs[2]
    s2 <- coefs[1]
    beta <- coefs[-(1:2)]
    SSE <- s2*n
    if (is.complex(eig)) 
        det <- sum(log(1 - rho * Re(eig)))
    else det <- sum(log(1 - rho * eig))
    ret <- (det - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
   ret
}

getVmat_eig <- function(coefs, n, eig, tol.solve=1.0e-10) {
    fd <- fdHess(coefs, f_laglm_eig, n, eig)
    solve(-(fd$Hessian), tol.solve=tol.solve)
}

f_laglm_Matrix <- function(coefs, n, W, I, nW, nChol, pChol) {
    rho <- coefs[2]
    beta <- coefs[-(1:2)]
    s2 <- coefs[1]
    SSE <- s2*n
    a <- -.Machine$double.eps^(1/2)
    b <- .Machine$double.eps^(1/2)
    Jacobian <- ifelse(rho > b, n * log(rho) +
        c(determinant(update(nChol, nW, 1/rho))$modulus),
        ifelse(rho < a, n* log(-(rho)) + c(determinant(update(pChol,
        W, 1/(-rho)))$modulus), 0.0))
    ret <- (Jacobian - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
   ret
}

getVmat_Matrix <- function(coefs, n, W, I, nW, nChol, pChol, tol.solve=1.0e-10) {
    fd <- fdHess(coefs, f_laglm_Matrix, n, W, I, nW, nChol, pChol)
    solve(-(fd$Hessian), tol.solve=tol.solve)
}

