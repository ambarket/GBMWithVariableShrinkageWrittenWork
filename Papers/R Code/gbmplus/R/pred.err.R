"pred.err" <-
function (lst) 
{
    dist <- lst[[1]]$distribution
    y <- lst[[1]]$mdata[,1]
    n <- length(y)
    if (dist == "gaussian") {
        pred.err <- mean((lst[[1]]$preds - y)^2)
        rel.pred.err <- pred.err/mean((y - mean(y))^2)
        list(pred.err = pred.err, rel.pred.err = rel.pred.err)
    }
    else if (dist == "laplace") {
        pred.err <- mean(abs(lst[[1]]$preds - y))
        rel.pred.err <- pred.err/mean(abs(y - mean(y)))
        list(pred.err = pred.err, rel.pred.err = rel.pred.err)
    }
    else if (dist == "bernoulli" || dist == "adaboost") {
        tab <- table(lst[[1]]$preds > 0, y)
        sdt <- sum(diag(tab))
        pred.err <- n - sdt
        rel.pred.err <- 1 - sdt/n
        list(pred.err = pred.err, rel.pred.err = rel.pred.err)
    }
    else cat(dist, " not yet implemented\n\n")
}

