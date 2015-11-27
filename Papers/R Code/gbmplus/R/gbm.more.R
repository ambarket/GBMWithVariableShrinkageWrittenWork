"gbm.more" <-
function(object,
                     n.new.trees = 100,
                     data = NULL,
                     weights = NULL,
                     offset = NULL,
                     verbose = NULL)
{
   if (!is.null(names(object))) {
        if(is.null(object$Terms) && is.null(object$data))
        {
            stop("The gbm model was fit using gbm.fit (rather than gbm) and keep.data was set to FALSE. gbm.more cannot locate the dataset.")
        }
        else if(is.null(object$data) && is.null(data))
        {
            stop("keep.data was set to FALSE on original gbm call and argument 'data' is NULL")
        }
        else if(is.null(object$data))
        {
            m <- eval(object$m, parent.frame())
            Terms <- attr(m, "terms")
            a <- attributes(Terms)
        
            y <- as.vector(model.extract(m, response))
            offset <- model.extract(m,offset)
            x <- model.frame(delete.response(Terms),
                            data,
                            na.action=na.pass)
        
            w <- weights
            if(length(w)==0) w <- rep(1, nrow(x))
            w <- w*length(w)/sum(w) # normalize to N
        
            if(is.null(offset) || (offset==0))
            {
                offset <- NA
            }
            Misc <- NA
        
            if(object$distribution == "coxph")
            {
                Misc <- as.numeric(y)[-(1:cRows)]
                y <- as.numeric(y)[1:cRows]
        
                # reverse sort the failure times to compute risk sets on the fly
                i.train <- order(-y[1:object$nTrain])
                i.test <- order(-y[(object$nTrain+1):cRows]) + object$nTrain
                i.timeorder <- c(i.train,i.test)
        
                y <- y[i.timeorder]
                Misc <- Misc[i.timeorder]
                x <- x[i.timeorder,,drop=FALSE]
                w <- w[i.timeorder]
                if(!is.na(offset)) offset <- offset[i.timeorder]
                object$fit <- object$fit[i.timeorder]
            }
        
            # create index upfront... subtract one for 0 based order
            if(ncol(x) > 1)
            {
                x.order <- apply(x[1:object$nTrain,],2,order,na.last=FALSE)-1
            }
            else
            {
                x.order <- order(x[1:object$nTrain,],na.last=FALSE)-1
            }
            x <- data.matrix(x)
            cRows <- nrow(x)
            cCols <- ncol(x)
        }
        else
        {
            y           <- object$data$y
            x           <- object$data$x
            x.order     <- object$data$x.order
            offset      <- object$data$offset
            Misc        <- object$data$Misc
            w           <- object$data$w
            cRows <- length(y)
            cCols <- length(x)/cRows
            if(object$distribution == "coxph")
            {
                i.timeorder <- object$data$i.timeorder
                object$fit <- object$fit[i.timeorder]
            }
        }
        
        if(is.null(verbose))
        {
            verbose <- object$verbose
        }
        x <- as.vector(x)
        
        gbm.obj <- .Call("gbm",
                        Y = as.double(y),
                        Offset = as.double(offset),
                        X = as.double(x),
                        X.order = as.integer(x.order),
                        weights = as.double(w),
                        Misc = as.double(Misc),
                        cRows = as.integer(cRows),
                        cCols = as.integer(cCols),
                        var.type = as.integer(object$var.type),
                        var.monotone = as.integer(object$var.monotone),
                        distribution = as.character(object$distribution),
                        n.trees = as.integer(n.new.trees),
                        interaction.depth = as.integer(object$interaction.depth),
                        n.minobsinnode = as.integer(object$n.minobsinnode),
                        shrinkage = as.double(object$shrinkage),
                        bag.fraction = as.double(object$bag.fraction),
                        nTrain = as.integer(object$nTrain),
                        fit.old = as.double(object$fit),
                        n.cat.splits.old = length(object$c.splits),
                        n.trees.old = as.integer(object$n.trees),
                        verbose = as.integer(verbose),
                        PACKAGE = "gbmplus")
        names(gbm.obj) <- c("initF","fit","train.error","valid.error",
                            "oobag.improve","trees","c.splits")
        
        gbm.obj$initF         <- object$initF
        gbm.obj$train.error   <- c(object$train.error, gbm.obj$train.error)
        gbm.obj$valid.error   <- c(object$valid.error, gbm.obj$valid.error)
        gbm.obj$oobag.improve <- c(object$oobag.improve, gbm.obj$oobag.improve)
        gbm.obj$trees         <- c(object$trees, gbm.obj$trees)
        gbm.obj$c.splits      <- c(object$c.splits, gbm.obj$c.splits)
        
        # cv.error not updated when using gbm.more
        gbm.obj$cv.error      <- object$cv.error
        
        gbm.obj$n.trees        <- length(gbm.obj$trees)
        gbm.obj$distribution   <- object$distribution
        gbm.obj$train.fraction <- object$train.fraction
        gbm.obj$shrinkage      <- object$shrinkage
        gbm.obj$bag.fraction   <- object$bag.fraction
        gbm.obj$var.type       <- object$var.type
        gbm.obj$var.monotone   <- object$var.monotone
        gbm.obj$var.names      <- object$var.names
        gbm.obj$interaction.depth <- object$interaction.depth
        gbm.obj$n.minobsinnode    <- object$n.minobsinnode
        gbm.obj$nTrain            <- object$nTrain
        gbm.obj$Terms             <- object$Terms
        gbm.obj$var.levels        <- object$var.levels
        gbm.obj$verbose           <- verbose
        
        gbm.obj$mdata <- object$m
        gbm.obj$nsamps <- object$nsamps
        gbm.obj$preds <- object$preds
        gbm.obj$type <- object$"cv"
        gbm.obj$m <- object$m.keep
        
        if(object$distribution == "coxph")
        {
            gbm.obj$fit[i.timeorder] <- gbm.obj$fit
        }
        if(!is.null(object$data))
        {
            gbm.obj$data <- object$data
        }
        else
        {
            gbm.obj$data <- NULL
            gbm.obj$m <- object$m
        }
        
   }     

   else {
        valid.err <- 0
        nsamps <- object[[1]]$nsamps
        n <- length(nsamps)
        gbm.obj <- list()
        
        for (i in 1:length(object)) {
          gbm.obj[[i]] <- gbm.more(object[[i]], n.new.trees = n.new.trees , 
          data = data, weights = weights, offset = offset, verbose = verbose)
          valid.err <- valid.err + gbm.obj[[i]]$valid.error*sum(nsamps==i)
        }

        valid.err <- valid.err/n
        n.trees <- gbm.obj[[1]]$n.trees
        valid.err <- smooth.spline(1:n.trees, valid.err)$y
        bb <- valid.err == min(valid.err)
        gbm.obj[[1]]$best <- best <- ((1:n.trees)[bb])[1]
        gbm.obj[[1]]$best.error <- best.error <- valid.err[best]
        gbm.obj[[1]]$mdata <- object[[1]]$mdata
        
        preds <- rep(NA,n)
        for (i in 1:length(object)) preds[nsamps==i] <- predict(gbm.obj[[i]], 
        newdata = gbm.obj[[1]]$mdata[nsamps==i, , drop = F], n.trees = best) 
        gbm.obj[[1]]$preds <- preds
    }
    class(gbm.obj) <- "gbm"
    return(gbm.obj)
   
}

