"summary.gbm" <-
function (object, plotit = TRUE, order = TRUE, n.trees = object$n.trees, 
    best = T, one = 1, colbar = "red", ...) 
{
    get.rel.inf <- function(obj) {
        lapply(split(obj[[6]], obj[[1]]), sum)
    }
    if (!is.null(names(object))) {
        vnames <- object$var.names
        if (object$n.trees < 1) {
            stop("n.trees must be greater than 0.")
        }
        if (n.trees > object$n.trees) {
            warning("Exceeded total number of GBM terms. Results use", 
                object$n.trees, "terms.\n")
            n.trees <- object$n.trees
        }
        temp <- unlist(lapply(object$trees[1:n.trees], get.rel.inf))
        rel.inf.compact <- unlist(lapply(split(temp, names(temp)), 
            sum))
        rel.inf.compact <- rel.inf.compact[names(rel.inf.compact) != 
            "-1"]
        rel.inf <- rep(0, length(vnames))
        i <- as.numeric(names(rel.inf.compact)) + 1
        rel.inf[i] <- rel.inf.compact
    }
    else {
        vnames <- object[[1]]$var.names
        nlist <- length(object)
        rel.inf.tot <- 0
        if (best) 
            best.n.tree <- object[[1]]$best
        else best.n.tree <- object[[1]]$n.trees
        for (i in 1:nlist) {
            temp <- unlist(lapply(object[[i]]$trees[one:best.n.tree], 
                get.rel.inf))
            rel.inf.compact <- unlist(lapply(split(temp, names(temp)), 
                sum))
            rel.inf.compact <- rel.inf.compact[names(rel.inf.compact) != 
                "-1"]
            rel.inf <- rep(0, length(vnames))
            i <- as.numeric(names(rel.inf.compact)) + 1
            rel.inf[i] <- rel.inf.compact
            rel.inf.tot <- rel.inf.tot + rel.inf
        }
        rel.inf <- rel.inf.tot/nlist
        rel.inf <- rel.inf.tot/nlist
    }
    if (order) {
        i <- order(-rel.inf)
    }
    else {
        i <- 1:length(rel.inf)
    }
    rel.inf <- 100 * rel.inf/sum(rel.inf)
    if (plotit) {
        barplot(rel.inf[rev(i)], horiz = TRUE, col = colbar, names = vnames[rev(i)], 
            xlab = "Relative influence", xpd = FALSE, ...)
    }
    return(data.frame(var = vnames[i], rel.influence = rel.inf[i]))
}

