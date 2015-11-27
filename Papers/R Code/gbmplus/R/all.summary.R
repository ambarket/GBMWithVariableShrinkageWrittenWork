"all.summary" <-
function (x, ..., pred = T, cor = T, xcor = T, sig = 3, samp = "auto") 
{
    if (!is.null(x$type) && x$type == "single") 
        stop("Single tree !!\n")
    else if (x[[1]]$type != "cv") 
        stop("Not a CV tree !!\n")
    old.par <- par(mfrow = c(2, 2), no.readonly = T)
    on.exit(par(old.par))
    xx <- x[[1]]
    res <- xx$mdata[, 1] - xx$preds
    plot(xx$preds, res, xlab = "Predicted", ylab = "Residuals")
    scatter.smooth(xx$preds, abs(res), xlab = "Predicted", ylab = "Absolute Residuals")
    qqnorm(res)
    cat("\nPredicted error\n")
    print(unlist(pred.err(x)), digits = sig + 1)
    cat("\n")
    sumry <- summary(x, ...)
    if (pred) {
        prd <- pred.vars(x)
        prd <- 100 * prd/sum(prd)
        ord <- match(sumry[, 1], names(prd))
        sumry <- cbind(sumry, prd[ord])
        colnames(sumry)[3] <- "partials"
    }
    if (cor) {
        cr <- corr.vars(x, samp = samp)
        if (!is.null(cr)) {
            ord <- match(sumry[, 1], names(cr[[1]]))
            sumry <- cbind(sumry, cr[[1]][ord])
            colnames(sumry)[3 + pred] <- "correlations"
            nr <- nrow(sumry)
            corpairs <- cr[[2]]
         }
    }
    sumry[, -1] <- zapsmall(sumry[, -1], digits = sig + 3)
    cat("Summary table\n")
    print(sumry, digits = sig)
    if (cor) {
        cat("\nPairwise correlations\n")
        print(corpairs, dig=3)
        cat("\nPairwise correlations cross-weighted by partials\n")
        print(sqrt(t(prd[-1]*t(prd[-nr] * corpairs^2))),dig=3)
    }
    cat("\n")
    if (xcor) {
        cat("Spearman correlations\n")
        cmat <- cor(x[[1]]$mdata, meth = "sp", use = "pairwise.complete.obs")
        cmat <- cmat[, -1]
        cmat <- cmat[-nrow(cmat), ]
        cmat[row(cmat) > col(cmat)] <- NA
        class(cmat) <- "anova"
        print(zapsmall(cmat, digits = sig))
        cat("\n")
    }
    if (cor) 
        invisible(list(sumry = sumry, corpairs = corpairs))
    else invisible(sumry)
}
