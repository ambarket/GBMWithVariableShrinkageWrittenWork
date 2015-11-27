"mvdsim" <-
function (n = 1000, nvar = 4, degree = 2, coefs = seq(1, 0.2, 
    length = degree), cor = 0, scale = T, err = 0.1, noise.vars = 0, b.scale = 1, 
    binary = FALSE) 
{
    x <- matrix(runif(n * nvar), ncol = nvar)
    if (cor) {
        for (i in 1:nvar) x[,i] <- x[,i] + cor * apply(x[,-1,drop=F]*runif(nvar-1,0.5,1),1,sum)
        x <- t(t(x) * sample(c(1,-1),nvar,rep=T))
    }
    x <- poly(x, degree = degree)
    deg <- attr(x, "degree")
    if (scale) 
        x <- scale(x)
    y <- x %*% (runif(length(deg), -1, 1) * coefs[deg])
    x <- x[, deg == 1]
    if (noise.vars) 
        x <- cbind(x, matrix(runif(noise.vars * n), nrow = n))
        x <- apply(x,2,function(x) (x-min(x))/(max(x)-min(x)))
    colnames(x) <- paste("x", 1:(nvar + noise.vars), sep = "")
    if (binary) {
        y <- y/sd(y) * b.scale
        ye <- rbinom(length(y), 1, 1/(1 + exp(-y)))
        data.frame(ye, y, x)
    }
    else {
        if (err) {
            err <- 1/(1/err - 1)
            y <- y/(sqrt(err) * sqrt(var(drop(y))))
            ye <- y + rnorm(n)
            data.frame(ye, y, x)
        }
        else data.frame(y, x)
    }
}
