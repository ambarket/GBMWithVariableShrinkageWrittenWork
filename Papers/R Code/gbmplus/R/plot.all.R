"plot.all" <-
function (obj, summary = TRUE, ord = T, n.vars, nn, main = NULL, ...) 
{
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    if (!summary) 
        ord <- FALSE
    if (!is.null(names(obj))) 
        vnames <- obj$var.names
    else vnames <- obj[[1]]$var.names
    if (missing(n.vars)) 
        n <- length(vnames)
    else n <- n.vars
    if (missing(nn)) 
        nn <- nbig.plt(n + summary)
    par(mfrow = nn, mar = c(5, 4, 1, 1) + 0.1)
    if (summary) 
        print(sumry <- summary.gbm(obj, main = main))
    if (ord) 
        ords <- match(sumry[, 1], vnames)[1:n]
    else ords <- 1:n
    for (i in ords) plot.gbm(obj, i, ...)
}
