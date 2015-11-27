"plot.gbm" <-
function (x, i.var = 1, npts, trim.grid = TRUE, number = 4, return.grid = FALSE, 
    center = TRUE, xlab, ylab, zlab = "Partial Dependence", type = c("link", 
        "response")[1], clip = TRUE, ptype = c("contour", "persp", 
        "lines", "trellis")[1], pcol = "vary", fpcol = heat.colors(50), 
    d = 6, theta = 40, n.trees, use = c("best", "all")[1], add.legend = TRUE, 
    se = FALSE, se.adjust = TRUE, se.smooth = c("mean", "smooth", 
        "none")[1], lty.se = 1, se.fill = TRUE, ylim, xrange = 0.025, 
    col1 = "red", col2 = "blue", pts = T, ...) 
{
    if (one.gbm <- !is.null(names(x))) {
        if (missing(n.trees)) {
            utype = c("best", "all")
            use <- utype[pmatch(use, utype)]
            if (use == "best") 
                n.trees <- x$best
            else n.trees <- x$n.trees
            if (is.null(n.trees)) 
                n.trees <- x$n.trees
        }
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
    if (!is.numeric(i.var)) {
        i.var <- unique(match(i.var, x$var.names))
    }
    if ((min(i.var) < 1) || (max(i.var) > length(x$var.names))) {
        warning("i.var must be between 1 and ", length(x$var.names))
    }
    if (length(i.var) > 3) {
        warning("plot.gbm creates up to 3-way interaction plots.\nplot.gbm will only return the plotting data structure.")
        return.grid = TRUE
    }
    se.sm <- c("mean", "smooth", "none")
    se.smooth <- se.sm[pmatch(se.smooth, c("mean", "smooth", 
        "none"))]
    pty <- c("contour", "persp", "lines", "trellis")
    ptype <- pty[pmatch(ptype, pty)]
    ty = c("link", "response")
    type <- ty[pmatch(type, ty)]
    distrib <- x$distribution
    if (missing(npts)) {
        if ((ptype == "contour" | ptype == "persp") & (length(i.var) < 
            3)) 
            npts <- rep(50, length(i.var))
        else npts <- c(50, rep(7, length(i.var) - 1))
    }
    if (length(trim.grid) == 1) 
        trim.grid <- rep(trim.grid, length(i.var))
    grid.levels <- vector("list", length(i.var))
    for (i in 1:length(i.var)) {
        if (is.numeric(x$var.levels[[i.var[i]]])) {
            grid.levels[[i]] <- seq(min(x$var.levels[[i.var[i]]]), 
                max(x$var.levels[[i.var[i]]]), length = npts[i])
            if (trim.grid[i]) 
                grid.levels[[i]] <- grid.levels[[i]][-c(1, npts[i])]
        }
        else {
            grid.levels[[i]] <- as.numeric(factor(x$var.levels[[i.var[i]]], 
                levels = x$var.levels[[i.var[i]]])) - 1
        }
    }
    X <- expand.grid(grid.levels)
    names(X) <- paste("x", 1:length(i.var), sep = "")
    facs <- all(!sapply(x$var.levels, is.character)[i.var][1:2])
    if (clip && (length(i.var) == 2) && facs) {
        xpdbox <- function(x, delta = 0.001) {
            x + delta * t(t(x) - apply(x, 2, mean))
        }
        poly.in <- function(xy, poly) {
            if (ncol(poly) == 2) 
                poly <- poly.tst(poly)
            n <- nrow(xy)
            np <- nrow(poly)
            nnp <- rep(n, np)
            check1 <- xor(xy[, 1] >= rep(poly[, 1], nnp), xy[, 
                1] > rep(poly[, 3], nnp))
            check2 <- rep(poly[, 2], nnp) + rep(poly[, 4], nnp) * 
                xy[, 1] > xy[, 2]
            as.vector(matrix(check1 & check2, n, np) %*% rep.int(1, 
                np)%%2 > 0)
        }
        poly.tst <- function(xy) {
            if (is.list(xy)) 
                xy <- as.data.frame(xy)
            if (!is.matrix(xy) || ncol(xy) != 2 || !is.numeric(xy[, 
                1]) || !is.numeric(xy[, 2])) 
                stop("xy must by nX2 numeric matrix or data.frame")
            xy <- as.matrix(xy)
            poly <- matrix(c(xy, xy[c(2:nrow(xy), 1), ]), ncol = 4)
            poly <- poly[poly[, 1] != poly[, 3], ]
            poly[, 4] <- (poly[, 4] - poly[, 2])/(poly[, 3] - 
                poly[, 1])
            poly[, 2] <- poly[, 2] - poly[, 1] * poly[, 4]
            poly
        }
        XX <- x$mdata[, x$var.names[i.var][1:2]]
        XX <- apply(XX, 2, function(x) {
            x[is.na(x)] <- mean(x, na.rm = T)
            x
        })
        hpts <- chull(XX)
        XXP <- xpdbox(XX[c(hpts, hpts[1]), ])
        ins <- poly.in(as.matrix(X), XXP)
    }
    if (one.gbm) {
        X$y <- .Call("gbm_plot", X = as.double(data.matrix(X)), 
            cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)), 
            i.var = as.integer(i.var - 1), n.trees = as.integer(n.trees), 
            initF = as.double(x$initF), trees = x$trees, c.splits = x$c.splits, 
            var.type = as.integer(x$var.type), PACKAGE = "gbmplus")
    }
    else if (!se) {
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
    else {
        ytemp <- matrix(NA, nrow = nrow(X), ncol = nlist)
        for (i in 1:nlist) {
            ytemp[, i] <- .Call("gbm_plot", X = as.double(data.matrix(X)), 
                cRows = as.integer(nrow(X)), cCols = as.integer(ncol(X)), 
                i.var = as.integer(i.var - 1), n.trees = as.integer(n.trees[i]), 
                initF = as.double(tree.list[[i]]$initF), trees = tree.list[[i]]$trees, 
                c.splits = tree.list[[i]]$c.splits, var.type = as.integer(x$var.type), 
                PACKAGE = "gbmplus")
        }
        X$y <- apply(ytemp, 1, mean, na.rm = T)
        X$se <- apply(ytemp, 1, var, na.rm = T)
        if (length(i.var) == 1 && is.numeric(x$var.levels[[i.var]])) {
            if (se.smooth == "smooth") 
                X$se <- runmed(X$se, k = 1 + 2 * floor(0.5 * 
                  sqrt(length(X$se))))
            else if (se.smooth == "mean") 
                X$se <- rep(mean(X$se), length(X$se))
            if (x$type == "cv") 
                X$se <- X$se * (nlist/2)
        }
        X$se <- sqrt(X$se)
    }
    if (clip && (length(i.var) == 2) && facs) 
        X$y[!ins] <- NA
    if (distrib == "bernoulli" && type == "response") 
        X$y <- 1/(1 + exp(-X$y))
    else if (distrib == "poisson" && type == "response") 
        X$y <- exp(X$y)
    else if (center) 
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
        if (se) 
            ylim <- c(min(X$y - se * X$se, na.rm = T), max(X$y + 
                se * X$se, na.rm = T))
        else ylim <- range(X$y, na.rm = T)
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
            plot(X$x1, X$y, type = "n", ylim = ylim, xlab = xlab, 
                ylab = ylab, ...)
            if (se) {
                if (!se.fill) {
                  lines(X$x1, X$y + se * X$se, col = col2, lty = lty.se)
                  lines(X$x1, X$y - se * X$se, col = col2, lty = lty.se)
                }
                else polygon(c(X$x1, rev(X$x1)), c(X$y + se * 
                  X$se, rev(X$y - se * X$se)), border = NA, col = "lightgray")
            }
            lines(X$x1, X$y, col = col1, ...)
            rug(x$var.levels[[i.var]])
        }
        else {
            if (!se) {
                barplot(X$y, names = x$var.levels[[i.var[i]]], 
                  ylim = ylim, xlab = xlab, ylab = ylab, col = col1, 
                  xpd = FALSE, ...)
                abline(h = 0)
            }
            else {
                yerr <- X$y + se * X$se * ifelse(X$y < 0, -1, 
                  1)
                z <- barplot(X$y, names = x$var.levels[[i.var[i]]], 
                  ylim = ylim, xlab = xlab, ylab = ylab, col = col1, 
                  xpd = FALSE, ...)
                abline(h = 0)
                segments(z, X$y, z, yerr)
                zerr <- 0.35 * (z[2] - z[1])
                segments(z + zerr, yerr, z - zerr, yerr)
            }
        }
    }
    else if (length(i.var) == 2) {
        if (!f.factor[1] && !f.factor[2]) {
            cp <- (ptype[1] == "persp") | (ptype[1] == "contour")
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[1]]
            if (missing(ylab) & cp) 
                ylab <- x$var.names[i.var[2]]
            if (missing(ylab) & !cp) 
                ylab <- paste("f(", x$var.names[i.var[1]], ",", 
                  x$var.names[i.var[2]], ")", sep = "")
            if (ptype[1] == "contour") {
                z <- matrix(X$y, ncol = length(grid.levels[[2]]))
                image(grid.levels[[1]], grid.levels[[2]], z, 
                  xlab = xlab, ylab = ylab, zlim = range(z, na.rm = T), 
                  xaxs = "r", yaxs = "r", col = fpcol, ...)
                contour(grid.levels[[1]], grid.levels[[2]], z, 
                  zlim = range(z, na.rm = T), add = T, col = col2, 
                  ...)
                if (pts) 
                  points(XX)
            }
            else if (ptype[1] == "persp") {
                z <- matrix(X$y, ncol = length(grid.levels[[2]]))
                if (pcol == "vary") 
                  pcol <- gcols(z, FUNCOL = fpcol)
                persp(grid.levels[[1]], grid.levels[[2]], z, 
                  zlab = zlab, xlab = xlab, ylab = ylab, zlim = range(z, 
                    na.rm = T), tick = "detailed", col = pcol, 
                  theta = theta, d = d, ...)
            }
            else if (ptype[1] == "lines") {
                X.lst <- split(X, X$x2)
                plot(X$x1, X$y, type = "n", xlab = xlab, ylab = ylab, 
                  ylim = ylim, ...)
                for (i in 1:length(X.lst)) lines(X.lst[[i]]$x1, 
                  X.lst[[i]]$y, col = i + 1)
                if (add.legend) 
                  legend(0.25 * min(X$x1) + 0.75 * max(X$x1), 
                    max(X$y, na.rm = T), as.character(signif(unique(sort(X$x2)), 
                      4)), lty = 1, col = 1 + (1:length(X.lst)), 
                    cex = 0.8, bty = "n")
            }
            else if (ptype[1] == "trellis") 
                if (!se) {
                  xyplot(y ~ x1 | x2, data = X, xlab = xlab, 
                    ylab = ylab, type = "l", ylim = ylim, ...)
                }
                else {
                  X$yhi <- X$y + se * X$se
                  X$ylo <- X$y - se * X$se
                  xyplot(y + yhi + ylo ~ x1 | x2, data = X, type = "l", 
                    allow.mult = T, xlab = xlab, ylab = ylab, 
                    col = c(2, 4, 4), lty = 1, ylim = ylim, ...)
                }
        }
        else if (f.factor[1] && !f.factor[2]) {
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[2]]
            if (missing(ylab)) 
                ylab <- paste("f(", x$var.names[i.var[2]], "|", 
                  x$var.names[i.var[1]], ")", sep = "")
            if (ptype[1] == "lines") {
                X.lst <- split(X, X$x1)
                plot(X$x2, X$y, type = "n", xlab = xlab, ylab = ylab, 
                  ylim = ylim, ...)
                for (i in 1:length(X.lst)) lines(X.lst[[i]]$x2, 
                  X.lst[[i]]$y, col = i + 1)
                if (add.legend) 
                  legend(0.25 * min(X$x2) + 0.75 * max(X$x2), 
                    max(X$y, na.rm = T), levels(X$x1), lty = 1, 
                    col = 1 + (1:length(X.lst)), cex = 0.8, bty = "n")
            }
            else if (!se) 
                xyplot(y ~ x1 | x2, data = X, xlab = xlab, ylab = ylab, 
                  type = "l", col = col1, ylim = ylim, ...)
            else {
                X$yhi <- X$y + se * X$se
                X$ylo <- X$y - se * X$se
                xyplot(y + yhi + ylo ~ x1 | x2, data = X, type = "l", 
                  allow.mult = T, xlab = xlab, ylab = ylab, col = c(col1, 
                    col2, col2), lty = 1, ylim = ylim, ...)
            }
        }
        else if (!f.factor[1] && f.factor[2]) {
            if (missing(xlab)) 
                xlab <- x$var.names[i.var[1]]
            if (missing(ylab)) 
                ylab <- paste("f(", x$var.names[i.var[1]], "|", 
                  x$var.names[i.var[2]], ")", sep = "")
            if (ptype[1] == "lines") {
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
            else if (!se) 
                xyplot(y ~ x1 | x2, data = X, xlab = xlab, ylab = ylab, 
                  type = "l", col = col1, ylim = ylim, ...)
            else {
                X$yhi <- X$y + se * X$se
                X$ylo <- X$y - se * X$se
                xyplot(y + yhi + ylo ~ x1 | x2, data = X, type = "l", 
                  allow.mult = T, xlab = xlab, ylab = ylab, col = c(col1, 
                    col2, col2), lty = 1, ylim = ylim, ...)
            }
        }
        else {
            if (ptype[1] == "trellis")
                barchart(x1 ~ y | x2, data = X, ylab = paste("f(~", 
                    x$var.names[i.var[1]], "*", x$var.names[i.var[2]], 
                    ")", sep = ""), orig = 0, ...)
            else {    
                if (!se) {
                    barplot(X$y, names = outer(substring(x$var.levels[[i.var[1]]],1,2),
                    substring(x$var.levels[[i.var[2]]],1,2),FUN=paste,sep="\n"),
                    ylim = ylim, xlab = xlab, ylab = ylab, col = col1, 
                    xpd = FALSE, ...)
                    abline(h = 0)
                }
                else {
                    yerr <- X$y + se * X$se * ifelse(X$y < 0, -1, 
                    1)
                    z <- barplot(X$y, names = outer(substring(x$var.levels[[i.var[1]]],1,2),
                    substring(x$var.levels[[i.var[2]]],1,2),FUN=paste,sep="\n"),
                    ylim = ylim, xlab = xlab, ylab = ylab, col = col1, 
                    xpd = FALSE, ...)
                    abline(h = 0)
                    segments(z, X$y, z, yerr)
                    zerr <- 0.35 * (z[2] - z[1])
                    segments(z + zerr, yerr, z - zerr, yerr)
                }
            }
        }
    }
    else if (length(i.var) == 3) {
        if (missing(xlab)) 
            xlab <- x$var.names[i.var[i[1]]]
        if (missing(ylab)) 
            ylab <- paste("f(", paste(x$var.names[i.var[1:3]], 
                collapse = ","), ")", sep = "")
        if (sum(f.factor) < 3) {
            if (!se) {
                xyplot(y ~ x1 | x2 * x3, data = X, xlab = xlab, 
                  ylab = ylab, type = "l", ylim = ylim, ...)
            }
            else {
                X$yhi <- X$y + se * X$se
                X$ylo <- X$y - se * X$se
                xyplot(y + yhi + ylo ~ x1 | x2 * x3, data = X, 
                  type = "l", allow.mult = T, xlab = xlab, col = c(2, 
                    4, 4), lty = 1, ylab = ylab, ylim = ylim, 
                  ...)
            }
        }
        else if (sum(f.factor) == 3) {
            barchart(x1 ~ y | x2 * x3, data = X, xlab = xlab, 
                ylab = ylab, origin = 0, ...)
        }
    }
}

