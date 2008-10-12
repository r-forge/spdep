aple.plot <- function(x, listw, do.plot=TRUE, ...) {
    pre <- spdep:::preAple(x=x, listw=listw)
    W2e <- eigen(pre$W2)
    SQRTW2 <- (diag(W2e$values^(0.5)) %*% t(W2e$vectors))
    X <- drop(SQRTW2 %*% x)
    NSQRTW2 <- (diag(W2e$values^(-0.5)) %*% t(W2e$vectors))
    Y <- drop(NSQRTW2 %*% pre$WU %*% x)
    plot(X, Y, ...)
    list(X=X, Y=Y)
}
