"gbm" <- function (formula = formula(data), data = list(), distribution = "gaussian", 
    weights, offset = NULL, var.monotone = NULL, n.trees = 200, 
    interaction.depth = 3, n.minobsinnode = ceiling(log10(nrow(data)) + 1),
    shrinkage = 0.05, bag.fraction = 0.5, keep.data = FALSE, verbose = FALSE,
    stratify = TRUE, monitor = TRUE, use = c("best","all")[1], seed = 0,
    na.omit.y = TRUE, cv.folds = 5, grp.col = 0, single.tree = FALSE, 
    train.fraction = 1) 
{
    m <- match.call(expand = FALSE)
    m$distribution <- m$weights <- m$offset <- m$var.monotone <- m$n.trees <- NULL
    m$interaction.depth <- m$n.minobsinnode <- m$shrinkage <- m$bag.fraction <- NULL
    m$keep.data <- m$verbose <- m$stratify <- m$monitor <- NULL
    m$use <- m$seed <- m$cv.folds <- m$grp.col <- m$single.tree <- m$train.fraction <- NULL
    m$na.omit.y <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m$na.action <- na.pass
    m.keep <- m
    m <- eval(m, parent.frame())
    Terms <- attr(m, "terms")
    a <- attributes(Terms)
    y <- model.extract(m, response)
    x <- model.frame(delete.response(Terms), data, na.action = na.pass)
    w <- model.extract(m, weights)
    if (na.omit.y & any(nays <- is.na(y))) {
    y <- y[!nays]
    x <- x[!nays,]
    w <- w[!nays]
    if (!is.null(offset)) offset <- offset[!nays]    
    }
    if (!is.null(model.extract(m, offset))) 
        stop("In gbm the offset term needs to be specified using the offset parameter.")
    var.names <- a$term.labels
    response.name <- dimnames(attr(terms(formula), "factors"))[[1]][1]
    nrows <- nrow(x)
    preds <- rep(NA, length = nrows)
    distributions <- c("bernoulli", "gaussian", "poisson", "adaboost", 
        "laplace", "coxph")
    distribution <- distributions[pmatch(distribution, distributions)]
    btype <- c("best", "all")
    use <- btype[pmatch(use, btype)]
    valid.err <- 0
    if (cv.folds > 0) {
        lst <- list()
        nontr <- list()
        if (seed > 0) 
            set.seed(seed)
        if (length(cv.folds) == nrows) {
            nsamps <- cv.folds
            cv.folds <- max(nsamps)
        }
        else if (grp.col == 0) {
            if (!stratify) 
                nsamps <- sample(rep(1:cv.folds, length = nrows))
            else nsamps <- rep(sample(1:cv.folds), ceiling(nrows/cv.folds))[1:nrows][order(y)]
        }
        else nsamps <- cv.var(data[, grp.col], cv.folds)
        for (i in 1:cv.folds) {
            nTrain <- length(trainees <- (1:nrows)[nsamps != 
                i])
            nontr[[i]] <- non.trainees <- (1:nrows)[is.na(match(1:nrows, 
                trainees))]
            temp.x <- x[c(trainees, non.trainees), , drop = F]
            temp.y <- y[c(trainees, non.trainees)]
            lst[[i]] <- gbm.fit(temp.x, temp.y, offset = offset, 
                distribution = distribution, w = w, var.monotone = var.monotone, 
                n.trees = n.trees, interaction.depth = interaction.depth, 
                n.minobsinnode = n.minobsinnode, shrinkage = shrinkage, 
                bag.fraction = bag.fraction, train.fraction = nTrain/nrows, 
                keep.data = keep.data, verbose = verbose, var.names = var.names, 
                response.name = response.name)
            valid.err <- valid.err + lst[[i]]$valid.error * length(non.trainees)
            lst[[i]]$Terms <- Terms
        }
        n.trees <- lst[[1]]$n.trees
        valid.err <- valid.err/nrows
        bb <- valid.err == min(valid.err)
        lst[[1]]$best <- best <- ((1:n.trees)[bb])[1]
        lst[[1]]$best.error <- best.error <- valid.err[best]
        for (i in 1:cv.folds) preds[nontr[[i]]] <- predict(lst[[i]], 
            newdata = x[nontr[[i]], , drop = F], n.trees = ifelse(use == 
                "best", best, n.trees))
        lst[[1]]$mdata <- m
        lst[[1]]$nsamps <- nsamps
        lst[[1]]$preds <- preds
        lst[[1]]$type <- "cv"
        lst[[1]]$m <- m.keep
    }
    if (single.tree || cv.folds == 0) {
        if (cv.folds > 0) 
            n.trees <- best
        lst <- gbm.fit(x = x, y = y, offset = offset, distribution = distribution, 
            w = w, var.monotone = var.monotone, n.trees = n.trees, 
            interaction.depth = interaction.depth, n.minobsinnode = n.minobsinnode, 
            shrinkage = shrinkage, bag.fraction = bag.fraction, 
            train.fraction = train.fraction, keep.data = keep.data, 
            verbose = verbose, var.names = var.names, response.name = response.name)
        lst$Terms <- Terms
        lst$best <- n.trees
        if (cv.folds > 0) 
            lst$best.error <- best.error
        lst$mdata <- m
        lst$type <- "single"
        lst$m <- m.keep
    }
    if (monitor) {
        cat("\nN trees =", n.trees, ": Depth =", interaction.depth, 
            ": Minimum node size =", n.minobsinnode, "\n")
        cat("Distribution =", distribution, ": Shrinkage =", 
            signif(shrinkage), ": Bag Fraction =", bag.fraction, 
            "\n")
        if (cv.folds) 
            cat("Best size =", best, ": Best error =", signif(best.error, 
                5), "\n\n")
        if (!single.tree & cv.folds) 
            monitor.gbm(lst)
    }
    class(lst) <- "gbm"
    return(lst)
}
