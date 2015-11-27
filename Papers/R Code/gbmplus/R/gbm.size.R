"gbm.size" <-
function(..., ints = c(1:5), FUN = gbm, auto.n.tree = 0) {
    old.par <- par(mfrow=c(2,1))
    on.exit(par(old.par))
    lints <- length(ints)
    trees <- list()
    if (auto.n.tree == 0)
        for (i in 1:lints) trees[[i]] <- FUN(..., int = i)
    else
        for (i in 1:lints) trees[[i]] <- FUN(..., int = i, n.tree = round(auto.n.tree/sqrt(i)))
    pred.err <- sapply(lapply(trees, pred.err),"[[",1)
    minpos <- (1:lints)[pred.err==min(pred.err)][1]
    plot(1:lints, pred.err, type="b", col="blue", xlab="Tree size", ylab="Predictied error", axes=F)
    axis(1,at=1:lints,lab=as.character(ints))
    axis(2)
    box()
    summary(trees[[minpos]])
    names(pred.err) <- as.character(ints)
    cat("Errors =",signif(pred.err,4),"\n")
    cat("Best depth of trees =",ints[minpos],": Minimum error =",signif(pred.err[minpos],4),"\n")
    trees[[minpos]]
    fit <- trees[[minpos]]
    fit$best.size <- ints[minpos]
    fit
}

