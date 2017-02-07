amat = matrix(c(1,0,1,0,0,1,0,1,1,0,1,1,0,1,1,1), nrow = 4, ncol = 4)
amat
data = replicate(100, rnorm(4))
S = matrix(0,4,4)
for (i in 1:100){
S = S + data[,i]%*%t(data[,i])
}
S = S/100;
rownames(S) <- c(1,2,3,4)
rownames(amat) <- c(1,2,3,4)
colnames(S) <- c(1,2,3,4)
colnames(amat) <- c(1,2,3,4)
nam <- rownames(S)
nod <- rownames(amat)
    if (is.null(nod)) {
        stop("The adjacency matrix has no labels.")
    }
if (!all(is.element(nod, nam))) 
    stop("The nodes of the graph do not match the names of the variables.")
else sek <- intersect(nam, nod)
S <- S[sek, sek, drop = FALSE]
amat <- amat[sek, sek, drop = FALSE]
nod <- rownames(amat)
if (all(amat == 0)) {
        alg <- 2
        cli = as.list(nod)
    }
cli = NULL
    if (is.null(cli)) {
        alg <- 3
    } else {
        alg <- 2
        nc <- length(cli)
        if (nc == 1) {
            return(list(Shat = S, dev = 0, df = 0, it = 1))}
    }
    k <- ncol(S)
alg
    if (alg == 1) {
        it <- 0
        W <- diag(diag(S))
        dimnames(W) <- dimnames(S)
        repeat {
            W.old <- W
            it <- it + 1
            for (i in 1:nc) {
                a <- cli[[i]]
                b <- setdiff(nod, a)
                Saa <- S[a, a]
                Waa <- W[a, a]
                Wba <- W[b, a]
                Wbb <- W[b, b]
                B <- Wba %*% solve(Waa)
                Spar <- Wbb - B %*% Waa %*% t(B)
                BV <- B %*% Saa
                W[b, a] <- BV
                W[a, b] <- t(BV)
                W[a, a] <- Saa
                W[b, b] <- Spar + BV %*% t(B)
            }
            if (sum(abs(W - W.old)) < tol) 
                break
        }
    }
    else if (alg == 2) {
        it = 0
        K <- solve(diag(diag(S)))
        dimnames(K) <- dimnames(S)
        repeat {
            K.old <- K
            it <- it + 1
            for (i in 1:nc) {
                a <- cli[[i]]
                b <- setdiff(nod, a)
                K[a, a] <- solve(S[a, a]) + K[a, b] %*% solve(K[b, 
                  b]) %*% K[b, a]
                if (pri) {
                  dev <- likGau(K, S, n, k)
                  cat(dev, "\n")
                }
            }
            if (sum(abs(K - K.old)) < tol) 
                break
        }
        W <- solve(K)
    }
pri = FALSE
    else if (alg == 3) {
        W0 <- S
        W <- S
        it <- 0
        converge = FALSE
        while (!converge) {
            it <- it + 1
            for (j in 1:k) {
                W11 <- W[-j, -j, drop = FALSE]
                w12 <- W[-j, j]
                s12 <- S[-j, j, drop = FALSE]
                paj <- amat[j, ] == 1
                paj <- paj[-j]
                beta <- rep(0, k - 1)
                if (all(!paj)) {
                  w <- rep(0, k - 1)
                } else {
                  beta[paj] <- solve(W11[paj, paj], s12[paj, 
                    ])
                  w <- W11 %*% beta}
                W[-j, j] <- w
                W[j, -j] <- w
            }
            di <- norm(W0 - W)
            if (pri) {
                cat(di, "\n")
            }
            if (di < tol) {
                converge <- TRUE
            } else {
                W0 <- W
            }
        }
    }
    df <- (sum(1 - amat) - k)/2
    Kh <- solve(W)
    dev <- likGau(Kh, S, n, k)
    list(Shat = W, dev = dev, df = df, it = it)
}
q()
