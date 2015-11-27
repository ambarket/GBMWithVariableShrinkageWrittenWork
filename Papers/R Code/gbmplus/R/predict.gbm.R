"predict.gbm" <-
function (object, newdata, n.trees, use = c("best","all")[1], type="link", 
         se = FALSE, se.adjust = TRUE, se.smooth = c("mean","smooth","none")[1], single.tree = FALSE, ...) 
{   
    if (missing(newdata)) stop("Specify newdata!!")
    ltype <- c("link","response")
    type <- ltype[pmatch(type, ltype)]
    ptype <- c("best","all")
    use <- ptype[pmatch(use, ptype)]
    one.gbm <- !is.null(names(object))
    if (!one.gbm) {
    objects <- object
    nlist <- length(objects)
    object <- objects[[1]]
    }
    if (!missing(n.trees)){
            change <- (n.trees>object$n.trees | n.trees<=0)
            if (any(change)) {
                n.trees[change] <- object$n.trees
                cat("Specified n.trees not valid: using n.trees =",n.trees,"\n")
            }
    }
    else if ((use == "best")&(!is.null(object$best))) n.trees <- object$best
    else n.trees <- object$n.trees 
    
    X <- model.frame(delete.response(object$Terms), newdata, 
        na.action = na.pass)
    cRows <- dim(X)[1]
    cCols <- dim(X)[2]
    for (i in 1:cCols) {
        if (is.factor(X[, i])) {
            j <- match(levels(X[, i]), object$var.levels[[i]])
            if (any(is.na(j))) {
                stop(paste("New levels for variable ", object$var.names[i], 
                  ": ", levels(X[, i])[is.na(j)], sep = ""))
            }
            X[, i] <- as.numeric(X[, i]) - 1
        }
    }
    X <- as.vector(unlist(X))
    if(missing(n.trees) || any(n.trees > object$n.trees))
    {
      n.trees[n.trees>object$n.trees] <- object$n.trees
      warning("Number of trees not specified or exceeded number fit so far. Using ",paste(n.trees,collapse=" "),".")
    }
    i.ntree.order <- order(n.trees)

    if (one.gbm) {
    predF <- .Call("gbm_pred", X = as.double(X), cRows = as.integer(cRows), 
        cCols = as.integer(cCols), n.trees = as.integer(n.trees), 
        initF = object$initF, trees = object$trees, c.split = object$c.split, 
        var.type = as.integer(object$var.type), single.tree = as.integer(single.tree),
        PACKAGE = "gbmplus")
    }    
    else  if (!se) {
    predF <- 0
    for (i in 1:nlist)
    predF <- predF + .Call("gbm_pred", X = as.double(X), cRows = as.integer(cRows), 
        cCols = as.integer(cCols), n.trees = as.integer(n.trees),
        initF = objects[[i]]$initF, trees = objects[[i]]$trees, c.split = objects[[i]]$c.split, 
        var.type = as.integer(objects[[i]]$var.type),  single.tree = as.integer(single.tree),
        PACKAGE = "gbmplus")   
    predF <- predF/nlist
    }
    else {
    ytemp <- matrix(NA, nrow = cRows, ncol = nlist)
    predF <- list()
    for (i in 1:nlist)
        ytemp[, i]  <- .Call("gbm_pred", X = as.double(X), cRows = as.integer(cRows), 
        cCols = as.integer(cCols), n.trees = as.integer(n.trees), 
        initF = objects[[i]]$initF, trees = objects[[i]]$trees, c.split = objects[[i]]$c.split, 
        var.type = as.integer(objects[[i]]$var.type),  single.tree = as.integer(single.tree),
        PACKAGE = "gbmplus")   
    predF$preds <- apply(ytemp, 1, mean, na.rm = T)
    predF$preds.se <- apply(ytemp, 1, var, na.rm = T) 
    if (se.smooth == "smooth")   
            predF$preds.se <- runmed(predF$preds.se, k=1+2*floor(0.5*sqrt(length(predF$se))))
    else if (se.smooth == "mean") 
            predF$preds.se <- rep(mean(predF$preds.se),length(predF$preds.se))
    if (object$type == "cv") 
        predF$preds.se <- predF$preds.se*(nlist/2)
    predF$preds.se <- sqrt(predF$preds.se)
    }


   if(type=="response")   {
      if(object$distribution=="bernoulli")  predF <- 1/(1+exp(-predF))
      else if(object$distribution=="poisson") predF <- exp(predF)
   }

   if(length(n.trees)>1)  {
      predF <- matrix(predF,ncol=length(n.trees),byrow=FALSE)
      colnames(predF) <- n.trees
      predF[,i.ntree.order] <- predF
   }

   if(!is.null(attr(object$Terms,"offset")))   {
      warning("predict.gbm does not add the offset to the predicted values.")
   }
   return(predF)
}

