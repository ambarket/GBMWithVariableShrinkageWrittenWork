"corr.vars" <-
function (x, samp = "auto", nn) 
{
    if (!is.null(x$type) && x$type == "single")
        stop("Single tree !!\n")
    else if (x[[1]]$type != "cv") 
        stop("Not a CV tree !!\n")
    n <- length(x[[1]]$fit)
    if (samp == "auto") 
        samp <- min(n, max(100, 3 * sqrt(n)))/n
    vars <- x[[1]]$var.names
    nvars <- length(vars)
    if (nvars == 1) 
        return()
    else if (x[[1]]$interaction.depth == 1) {
        res <- rep(0, nvars)
        names(res) <- vars
        resmat <- matrix(0, nrow = nvars - 1, ncol = nvars - 
            1)
        resmat[row(resmat) > col(resmat)] <- NA
        dimnames(resmat) <- list(vars[-nvars], vars[-1])
        class(resmat) <- "anova"
        list(res, resmat)
    }
    else {
        data <- x[[1]]$mdata[, vars]
        if (samp & samp < 1) {
            cat("Correlations subsampling fraction =", signif(samp, 
                3), "\n")
            samp <- round(samp * n)
            data <- data[sample(n, samp), ]
        }
        nsamps <- x[[1]]$nsamps
        nnsamps <- max(nsamps)
        if (missing(nn)) {
            res <- rep(NA, nvars)
            names(res) <- vars
            resmat <- matrix(NA, nrow = nvars - 1, ncol = nvars - 
                1)
            dimnames(resmat) <- list(vars[-nvars], vars[-1])
            pr.ones <- pr.comp <- list()
            pr.all <- peffects.gbm(x, data = data, ret = T, cent = F)[, 
                nvars + 1]
            for (i in 1:nvars) {
                pr.comp[[i]] <- peffects.gbm(x, data = data[, 
                  (1:nvars)[-i], drop = F], ret = T, cent = F)[, 
                  nvars]
                pr.ones[[i]] <- peffects.gbm(x, data = data[, 
                  (1:nvars)[i], drop = F], ret = T, cent = F)[, 
                  2]
                res[i] <- sqrt(1 - summary(lm(pr.all ~ pr.ones[[i]] + 
                  pr.comp[[i]]))$r.sq)
            }
            for (i in 1:(nvars - 1)) {
                for (j in (i + 1):nvars) {
                  pr.pair <- peffects.gbm(x, data = data[, c(i, 
                    j)], ret = T, cent = F)[, 3]
                  resmat[i, j - 1] <- sqrt(1 - summary(lm(pr.pair ~ 
                    pr.ones[[i]] + pr.ones[[j]]))$r.sq)
                }
                class(resmat) <- "anova"
            }
            list(res, resmat)
        }
        else {
            nx <- length(nn)
            preds1 <- peffects.gbm(x, data = data, ret = T, cent = F)[, 
                nvars + 1]
            preds2 <- peffects.gbm(x, data = data[, (1:nvars)[-nn], 
                drop = F], ret = T, cent = F)[, nvars - nx + 
                1]
            preds3 <- peffects.gbm(x, data = data[, (1:nvars)[nn], 
                drop = F], ret = T, cent = F)[, nx + 1]
            res <- sqrt(1 - summary(lm(preds1 ~ preds2 * preds3))$r.sq)
            res
        }
    }
}
