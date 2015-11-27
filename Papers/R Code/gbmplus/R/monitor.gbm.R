"monitor.gbm" <-
function (obj) 
{
    n<- length(obj)
    nt <- 1:obj[[1]]$n.trees
    ind <- unique(round(nt, -floor(log(nt, 10))))
    ind[length(ind)] <- max(nt)
    t1 <- t2 <- matrix(NA,nrow=length(ind),ncol=n)
    for (i in 1:n) {
    t1[,i] <- obj[[i]]$train.error[ind] 
    t2[,i] <- obj[[i]]$valid.error[ind] 
    }
    tr.err <- apply(t1,1, mean)
    val.err <- apply(t2,1, mean)
    df <- cbind(TrainLL = tr.err, ValidLL = val.err)
    row.names(df) <- ind
    print(df)
    cat("\n")
}

