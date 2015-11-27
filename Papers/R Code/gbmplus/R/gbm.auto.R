"gbm.auto" <-
function (form, data, ..., pred = 1, cor.keep = 0.25, mono = F, 
    size = F, cors = T, seed = sample(1:100, 1), monitor = F, 
    np = c(5, 3)[1], long = F) 
{
    fit <- gbm(form, data, ...)
    pred.err <- pred.err(fit)
    cat("\nPredicted error\n")
    print(unlist(pred.err))
    pred.vars <- pred.vars(fit)
    cat("\nPredictors \n")
    print(pred.vars)
    pred <- pred * sum(pred.vars) * pred.err[[2]]/length(fit[[1]]$var.names)
    cat("\nCritical value for dropping =", signif(pred, 3), "\n")
    keeps <- pred.vars > pred
    if (cors) {
        cor.vars <- corr.vars(fit)
        cat("\nCorrelations of one vs remaining predictors\n")
        print(cor.vars[[1]])
        cat("\nPairwise correlations\n")
        print(cor.vars[[2]])
        cat("\nCritical value for dropping =", signif(cor.keep, 
            3), "\n")
        keeps <- keeps | (cor.vars[[1]] > cor.keep)
    }
    xnames <- fit[[1]]$var.names[keeps]
    if (all(keeps)) 
        cat("\nNo predictors dropped\n")
    else {
        cat("\nDropping: ", fit[[1]]$var.names[!keeps], "\n")
        dropnames <- paste(fit[[1]]$var.names[!keeps], collapse = " - ")
        newform <- formula(paste("~.-", dropnames, collapse = ""))
        form <- update.formula(form, newform)
    }
    print(form)
    fit <- gbm(form, data, monitor=F, seed = seed, 
        ...)
    pred.err <- pred.err(fit)
    cat("\nPredicted error\n")
    print(unlist(pred.err))
    pred.vars <- pred.vars(fit)
    cat("\nPredictors \n")
    print(pred.vars)
    nxs <- length(xnames)
    if (nxs > 1) {
        best.pred.err <- pred.err[[2]]
        stop.now <- FALSE
        best <- 0
        xnames <- fit[[1]]$var.names
        while (!stop.now && length(xnames)>1) {
            cat("\nCurrent pedictors:", xnames, "\n")
            perrs <- list()
            for (i in 1:nxs) {
                if (long) 
                  cat("\nTest drop", xnames[i], "\n")
                tempform <- formula(paste("~.-", xnames[i], collapse = ""))
                dropform <- update.formula(form, tempform)
                if (long) 
                  print(dropform)
                fit <- gbm(dropform, data, monitor=F, seed = seed, ...)
                perrs[[i]] <- pred.err(fit)
                if (long) {
                  cat("Predicted error\n")
                  print(unlist(perrs[[i]]))
                }
            }
            perrs1 <- sapply(perrs, "[[", 2)
            best <- min(perrs1)
            nbest <- ((1:nxs)[perrs1 == best])[1]
            if (best < best.pred.err) {
                best.pred.err <- best
                stop.now <- FALSE
                newform <- formula(paste("~.-", xnames[nbest], 
                  collapse = ""))
                form <- update.formula(form, newform)
                xnames <- xnames[-nbest]
                nxs <- length(xnames)
            }
            else stop.now <- TRUE
        }
    }
    if (size) {
        cat("\nEvaluating best size\n\n")
        best.size <- gbm.size(form, data, monitor=F, seed = seed, 
            ...)$best.size
        fit <- gbm(form, data, int = best.size, monitor = F, 
            seed = seed, ...)
    }
    if (size) 
        fit <- gbm(form, data, int = best.size, monitor = F, 
            seed = seed, ...)
    else fit <- gbm(form, data, monitor = F, seed = seed, 
        ...)
    pred.err <- pred.err(fit)
    cat("\nPredicted error\n")
    print(unlist(pred.err))
    if (mono) {
        if (np == 5) 
            levs <- c(2, 4, 6, 8, 10)
        else if (np == 3) 
            levs <- c(2, 6, 10)
        vtype <- fit[[1]]$var.type
        nxs <- length(vtype)
        vmono <- rep(0, nxs)
        for (i in 1:nxs) {
            if (vtype[i] == 0) {
                xdata <- data.frame(fit[[1]]$var.levels[[i]][levs])
                names(xdata) <- fit[[1]]$var.names[i]
                ord <- order(peffects.gbm(fit, data = xdata, ret = T)[, 2])
                if (all(ord == (1:np))) vmono[i] <- 1
                else if (all(ord == (np:1))) vmono[i] <- -1
            }
        }
        cat("\nMonotone predictors", vmono, "\n")
        if (any(vmono != 0)) {
            if (size) 
                fit <- gbm(form, data, int = best.size, monitor=F, 
                  var.mono = vmono, seed = seed, ...)
            else fit <- gbm(form, data, var.mono = vmono, monitor=F,
                seed = seed, ...)
        }
    }
    pred.vars <- pred.vars(fit)
    cat("\nPredictors \n")
    print(pred.vars)
    all.summary(fit, plot = T, cor = cors)
    cat("Final model:\n")
    print(form)
    cat("\n")
    fit
}

