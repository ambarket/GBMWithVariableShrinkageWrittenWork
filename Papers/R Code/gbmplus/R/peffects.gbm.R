"peffects.gbm" <-
function (x, i.var = 1, data, trim.grid = TRUE, number = 4, return.grid = FALSE, 
    center = TRUE, ngrp = 4, xlab, ylab, zlab = "Partial Dependence", 
    stype = c("logit", "probability"), clip = TRUE, ptype = c("contour", 
        "persp", "lines", "trellis")[1], pcols = "white", n.trees, 
    use = c("best", "all")[1], add.legend = TRUE, ylim, xrange = 0.025, 
    col1 = 1, col2 = 2, pts = T, ...) 
{
    if (one.gbm <- !is.null(names(x))) {
        if (missing(n.trees)) {
            utype = c("best", "all")
            use <- utype[pmatch(use, utype)]
            if (use == "best") 
                n.trees <- x$best
            else n.trees <- x$n.trees
        }
        else n.trees <- x$n.trees
        se <- FALSE
    }
    else {
        tree.list <- x
        nlist <- length(tree.list)
        x <- tree.list[[1]]
        if (missing(n.trees)) {
            utype = c("best", "all")
            use <- utype[pmatch(use, utype)]
            if (use == "best") 
                n.trees <- rep(x$best, nlist)
            else n.trees <- rep(x$n.trees, nlist)
        }
        else n.trees <- rep(x$n.trees, nlist)
    }
    if (missing(data)) {
        if (!is.numeric(i.var)) {
            i.var <- unique(match(i.var, x$var.names))
        }
        if ((min(i.var) < 1) || (max(i.var) > length(x$var.names))) {
            warning("i.var must be between 1 and ", length(x$var.names))
        }
        if (length(i.var) > 3) {
            return.grid = TRUE
        }
        data <- x$mdata[x$var.names]
        nxvar <- length(x$var.names)
        nrows <- nrow(data)
        X <- data[, i.var, drop = F]
    }
    else {
        i.var <- (1:length(x$var.names))[x$var.names %in% names(data)]
        nrows <- nrow(data)
        X <- data
    }
    for (i in 1:ncol(X)) {
        if (is.factor(X[, i])) 
            X[, i] <- number(X[, i]) - 1
    }
    if (length(i.var) > 1) 
        colnames(X) <- paste("x", 1:length(i.var), sep = "")
    else colnames(X) <- "x1"
    pty <- c("contour", "persp", "lines", "trellis")
    ptype <- pty[pmatch(ptype, pty)]
    distrib <- x$distribution
    lcols <- trellis.par.get()$superpose.line$col
    col1 <- lcols[col1]
    col2 <- lcols[col2]
    if (one.gbm) {
        X$y <- .Call("gbm_plot", X = as.double(data.matrix(X)), 
            cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)), 
            i.var = as.integer(i.var - 1), n.trees = as.integer(n.trees), 
            initF = as.double(x$initF), trees = x$trees, c.splits = x$c.splits, 
            var.type = as.integer(x$var.type), PACKAGE = "gbmplus")
    }
    else {
        ytemp <- 0
        for (i in 1:nlist) {
            ytemp <- ytemp + .Call("gbm_plot", X = as.double(data.matrix(X)), 
                cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)), 
                i.var = as.integer(i.var - 1), n.trees = as.integer(n.trees[i]), 
                initF = as.double(tree.list[[i]]$initF), trees = tree.list[[i]]$trees, 
                c.splits = tree.list[[i]]$c.splits, var.type = as.integer(x$var.type), 
                PACKAGE = "gbmplus")
        }
        X$y <- ytemp/nlist
    }
    if (center) 
        X$y <- X$y - mean(X$y, na.rm = T)
    f.factor <- rep(FALSE, length(i.var))
    for (i in 1:length(i.var)) {
        if (!is.numeric(x$var.levels[[i.var[i]]])) {
            X[, i] <- factor(x$var.levels[[i.var[i]]][X[, i] + 
                1], levels = x$var.levels[[i.var[i]]])
            f.factor[i] <- TRUE
        }
    }
    if (return.grid) {
        names(X)[1:length(i.var)] <- x$var.names[i.var]
        return(X)
    }
    nrows <- nrow(x$mdata)
    if (missing(ylim)) {
        ylim <- range(X$y, na.rm = T)
        if (xrange) {
            rnge <- diff(ylim)
            ylim[1] <- ylim[1] - xrange * rnge
            ylim[2] <- ylim[2] + xrange * rnge
        }
    }
    if (length(i.var) == 1) {
        if (missing(xlab)) 
            xlab <- x$var.names[i.var]
        if (missing(ylab)) 
            ylab <- paste("f(", x$var.names[i.var], ")", sep = "")
        if (!f.factor) {
            X <- X[order(X$x1), ]
            plot(X$x1, X$y, type = "n", ylim = ylim, xlab = xlab, 
                ylab = ylab, ...)
            lines(X$x1, X$y, col = col1, ...)
            rug(x$var.levels[[i.var]])
        }
        else {
            X <- X[!duplicated(X$x1), ]
            X <- X[order(X$x1), ]
            barplot(X$y, names = x$var.levels[[i.var[i]]], ylim = ylim, 
                xlab = xlab, ylab = ylab, col = col1, xpd = FALSE, 
                ...)
            abline(h = 0)
        }
    }
    else if (length(i.var) == 2) {
        if (!f.factor[1] && !f.factor[2]) {
            cp <- (ptype[1] == "persp") | (ptype[1] == "contour")
            if (cp) require(akima)
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[1]]
            if (missing(ylab) & cp) 
                ylab <- x$var.names[i.var[2]]
            if (missing(ylab) & !cp) 
                ylab <- paste("f(", x$var.names[i.var[1]], ",", 
                  x$var.names[i.var[2]], ")", sep = "")
            if (ptype[1] == "contour") {
                z <- interp(X$x1, X$x2, X$y, duplicate = "mean")
                image(z, xlab = x$var.names[i.var[1]], ylab = x$var.names[i.var[2]], 
                  zlim = range(z$z, na.rm = T), col = heat.colors(100))
                contour(z, zlim = range(z$z, na.rm = T), add = T, 
                  col = col2)
                if (pts) 
                  points(X$x1, X$x2)
            }
            else if (ptype[1] == "persp") {
                z <- interp(X$x1, X$x2, X$y, duplicate = "mean")
                persp(z, zlab = "partial dependence", xlab = x$var.names[i.var[1]], 
                  ylab = x$var.names[i.var[2]], zlim = range(z$z, 
                    na.rm = T), col = pcols, tick = "detailed", 
                  ...)
            }
            else if (ptype[1] == "lines") {
                X <- X[order(X$x1), ]
                X$x2 <- cut(X$x2, breaks = quantile(X$x2, prob = (0:ngrp)/ngrp), 
                  inc = T)
                X.lst <- split(X, X$x2)
                plot(X$x1, X$y, type = "n", xlab = xlab, ylab = ylab, 
                  ylim = ylim, ...)
                for (i in 1:length(X.lst)) lines(X.lst[[i]]$x1, 
                  X.lst[[i]]$y, col = i + 1)
                if (add.legend) 
                  legend(0.25 * min(X$x1) + 0.75 * max(X$x1), 
                    max(X$y, na.rm = T), levels(X$x2), lty = 1, 
                    col = 1 + (1:length(X.lst)), cex = 0.8, bty = "n")
            }
            else if (ptype[1] == "trellis") {
                X <- X[order(X$x1), ]
                X$x2 <- equal.count(X$x2, ngrp, overlap = 0)
                xyplot(y ~ x1 | x2, data = X, xlab = xlab, ylab = ylab, 
                  type = "l", ylim = ylim, ...)
            }
        }
        else if (f.factor[1] && !f.factor[2]) {
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[2]]
            if (missing(ylab)) 
                ylab <- paste("f(", x$var.names[i.var[2]], "|", 
                  x$var.names[i.var[1]], ")", sep = "")
            X <- X[order(X$x2), ]
            X.lst <- split(X, X$x1)
            if (ptype[1] == "lines") {
                plot(X$x2, X$y, type = "n", xlab = xlab, ylab = ylab, 
                  ylim = ylim, ...)
                for (i in 1:length(X.lst)) lines(X.lst[[i]]$x2, 
                  X.lst[[i]]$y, col = i + 1)
                if (add.legend) 
                  legend(0.25 * min(X$x2) + 0.75 * max(X$x2), 
                    max(X$y, na.rm = T), levels(X$x1), lty = 1, 
                    col = 1 + (1:length(X.lst)), cex = 0.8, bty = "n")
            }
            else xyplot(y ~ x2 | x1, data = X, xlab = xlab, ylab = ylab, 
                type = "l", col = col1, ylim = ylim, ...)
        }
        else if (!f.factor[1] && f.factor[2]) {
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[1]]
            if (missing(ylab)) 
                ylab <- paste("f(", x$var.names[i.var[1]], "|", 
                  x$var.names[i.var[2]], ")", sep = "")
            X <- X[order(X$x1), ]
            X.lst <- split(X, X$x2)
            if (ptype[1] == "lines") {
                plot(X$x1, X$y, type = "n", xlab = xlab, ylab = ylab, 
                  ylim = ylim, ...)
                for (i in 1:length(X.lst)) lines(X.lst[[i]]$x1, 
                  X.lst[[i]]$y, col = i + 1)
                if (add.legend) 
                  legend(0.25 * min(X$x1) + 0.75 * max(X$x1), 
                    max(X$y, na.rm = T), levels(X$x2), lty = 1, 
                    col = 1 + (1:length(X.lst)), cex = 0.8, bty = "n")
            }
            else xyplot(y ~ x1 | x2, data = X, xlab = xlab, ylab = ylab, 
                type = "l", col = col1, ylim = ylim, ...)
        }
        else {
            barchart(x1 ~ y | x2, data = X, ylab = paste("f(~", 
                x$var.names[i.var[1]], "*", x$var.names[i.var[2]], 
                ")", sep = ""), orig = 0, ...)
        }
    }
    else if (length(i.var) == 3) {
        if (missing(xlab)) 
            xlab <- x$var.names[i.var[i[1]]]
        if (missing(ylab)) 
            ylab <- paste("f(", paste(x$var.names[i.var[1:3]], 
                collapse = ","), ")", sep = "")
        if (sum(f.factor) < 3) {
            xyplot(y ~ x1 | x2 * x3, data = X, xlab = xlab, ylab = ylab, 
                type = "l", ylim = ylim, ...)
        }
        else if (sum(f.factor) == 3) {
            barchart(x1 ~ y | x2 * x3, data = X, xlab = xlab, 
                ylab = ylab, ...)
        }
    }
}

