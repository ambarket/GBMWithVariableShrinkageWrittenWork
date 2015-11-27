"pred.vars" <-
function(x) {

#   Returns a measure of how well each x-var predicts
#   The value is for x1 : SS(pred(x1)-mean(y)) / SS(pred(allx)-meany)
#   Xvars with values close to zero (<0.001 ???) can be dropped ?????
#
    if (!is.null(x$type) && x$type == "single")
        stop("Single tree !!\n")
    else  if (x[[1]]$type!="cv") 
        stop("Not a CV tree !!\n")
    mse <- mean((x[[1]]$preds - mean(x[[1]]$mdata[,1]))^2)
    vars <- x[[1]]$var.names
    nvars <- length(vars)
    data <- x[[1]]$mdata[, -1, drop = F]
    n <- nrow(data)
    preds <- matrix(NA, nrow = n, ncol = nvars)
    nsamps <- x[[1]]$nsamps
    nnsamps <- max(nsamps)
    for (i in 1:nvars) for (j in 1:nnsamps) preds[nsamps == j, 
        i] <- peffects.gbm(x[[j]], data = data[nsamps == j, i, 
        drop = F], ret = T, cent = F, n.tree = x[[1]]$best)[, 2]
    preds <- preds - mean(x[[1]]$mdata[,1])
    res <- apply(preds^2, 2, mean)/mse
    names(res) <- vars
    res 
}

