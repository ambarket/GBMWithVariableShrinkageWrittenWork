"plot.err" <-
function (obj, add.lines = T) 
{
    if (class(obj) != "gbm") 
        stop("Not GBM object\n")
    if (!is.null(names(obj))) {
        ylims <- range(c(obj$valid.error, obj$train.error))
        plot(obj$train.error, type = "l", xlab = "N Trees", ylab = "Train/Valid Error")
        lines(obj$valid.error, type = "l", col = 2)
        best.tree <- (1:obj$n.tree)[obj$valid.error == min(obj$valid.error)][1]
        best.valid.error <- obj$valid.error[best.tree]
        cat("Best size tree =", best.tree, " : Best valid error =", 
            signif(best.valid.error, 4), "\n")
    }
    else {
        nlist <- length(obj)
        n <- length(obj[[1]]$valid.error)
        best <- obj[[1]]$best
        best.error <- obj[[1]]$best.error
        mat <- matrix(NA, ncol = nlist + 1, nrow = n)
        best.ns <- rep(NA, nlist)
        bests <- rep(NA, nlist)
        distrib <- obj[[1]]$distribution
        if ((distrib == "gaussian") || (distrib == "adaboost") || 
            (distrib == "laplace") || (distrib == "bernoulli")) 
            FUN <- min
        else FUN <- max
        for (i in 1:nlist) {
            mat[, i + 1] <- obj[[i]]$valid.error
            best.ns[i] <- ((1:n)[mat[, i + 1] == FUN(mat[, i + 
                1])])[1]
            bests[i] <- mat[best.ns[i], i + 1]
        }
        mat[, 1] <- apply(mat[, -1], 1, mean)
        
        matplot(1:n, mat, col = 1:(nlist + 1), type = "l", lty = 1, 
            lwd = c(2, rep(1, nlist)), xlab = "N trees", ylab = "Mean valid error")
        if (add.lines) {
            abline(v = best, lty = 2, col = "gray")
            abline(h = best.error, lty = 2, col = "gray")
        }
        cat("Best ave min n's =", best.ns, "\n")
        cat("Best ave pred errs =", signif(bests, 3), "\n")
        cat("Best min ave pred err =", signif(best.error, 3), " : best ave min n tree =", 
            best, "\n\n")
        invisible(list(best.error = best.error, best = best))
    }
}

