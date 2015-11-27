"gbm.fit" <-
function(x,y,
                    offset = NULL,
                    misc = NULL,
                    distribution = "gaussian",
                    w = NULL,
                    var.monotone = NULL,
                    n.trees = 200,
                    interaction.depth = 3,
                    n.minobsinnode = ceiling(log10(nrow(data))+1),
                    shrinkage = 0.1,
                    bag.fraction = 0.5,
                    train.fraction = 1.0,
                    keep.data = TRUE,
                    verbose = TRUE,
                    var.names = NULL,
                    response.name = NULL)
{
   cRows <- nrow(x)
   cCols <- ncol(x)
   if(is.null(var.names))
   {
      if(is.matrix(x))
      {
         var.names <- colnames(x)
      }
      else if(is.data.frame(x))
      {
         var.names <- names(x)
      }
      else
      {
         var.names <- paste("X",1:cCols,sep="")
      }
   }
   if(is.null(response.name))
   {
      response.name <- "y"
   }

   # check dataset size
   if(cRows*train.fraction*bag.fraction <= 2*n.minobsinnode+1)
   {
      stop("The dataset size is too small or subsampling rate is too large: cRows*train.fraction*bag.fraction <= n.minobsinnode")
   }

   if(nrow(x) != ifelse(class(y)=="Surv", nrow(y), length(y)))
   {
      stop("The number of rows in x does not equal the length of y.")
   }

   if(interaction.depth < 1)
   {
      stop("interaction.depth must be at least 1.")
   }

   if(length(w)==0) w <- rep(1, cRows)
   else if(any(w < 0)) stop("negative weights not allowed")
   w <- w*length(w)/sum(w) # normalize to N

   if(any(is.na(y))) stop("Missing values are not allowed in the response, ",
                          response.name)

   if(is.null(offset) || (offset==0))
   {
      offset <- NA
   }
   else if(length(offset) != length(y))
   {
      stop("The length of offset does not equal the length of y.")
   }
   Misc <- NA

   # setup variable types
   var.type <- rep(0,cCols)
   var.levels <- vector("list",cCols)
   for(i in 1:length(var.type))
   {
      if(all(is.na(x[,i])))
      {
         stop("variable ",i,": ",var.names[i]," has only missing values.")
      }
      if(is.ordered(x[,i]))
      {
         var.levels[[i]] <- levels(x[,i])
         x[,i] <- as.numeric(x[,i])-1
         var.type[i] <- 0
      }
      else if(is.factor(x[,i]))
      {
         if(length(levels(x[,i]))>1024)
            stop("gbm does not currently handle categorical variables with more than 1024 levels. Variable ",i,": ",var.names[i]," has ",length(levels(x[,i]))," levels.")
         var.levels[[i]] <- levels(x[,i])
         x[,i] <- as.numeric(x[,i])-1
         var.type[i] <- max(x[,i],na.rm=TRUE)+1
      }
      else if(is.numeric(x[,i]))
      {
         var.levels[[i]] <- quantile(x[,i],prob=(0:10)/10,na.rm=TRUE)
      }
      else
      {
         stop("variable ",i,": ",var.names[i]," is not of type numeric, ordered, or factor.")
      }

      # check for some variation in each variable
      if(length(unique(var.levels[[i]])) == 1)
      {
         warning("variable ",i,": ",var.names[i]," has no variation.")
      }
   }

   nTrain <- as.integer(train.fraction*cRows)

   supported.distributions <-
      c("bernoulli","gaussian","poisson","adaboost","laplace","coxph")
   # check potential problems with the distributions
   if(!is.element(distribution,supported.distributions))
   {
      stop("Distribution ",distribution," is not supported")
   }
   if((distribution == "bernoulli") && !all(is.element(y,0:1)))
   {
      stop("Bernoulli requires the response to be in {0,1}")
   }
   if((distribution == "poisson") && any(y<0))
   {
      stop("Poisson requires the response to be positive")
   }
   if((distribution == "poisson") && any(y != trunc(y)))
   {
      stop("Poisson requires the response to be a positive integer")
   }
   if((distribution == "adaboost") && !all(is.element(y,0:1)))
   {
      stop("This version of AdaBoost requires the response to be in {0,1}")
   }
   if((distribution == "laplace") && (length(unique(w)) > 1))
   {
      stop("This version of gbm for the Laplace loss lacks a weighted median. For now the weights must be constant.")
   }
   if(distribution == "coxph")
   {  
      require(survival)  
      if(class(y)!="Surv")
      {
         stop("Outcome must be a survival object Surv(time,failure)")
      }
      if(attr(y,"type")!="right")
      {
         stop("gbm() currently only handles right censored observations")
      }
      Misc <- y[,2]
      y <- y[,1]

      # reverse sort the failure times to compute risk sets on the fly
      i.train <- order(-y[1:nTrain])
      n.test <- cRows - nTrain
      if(n.test > 0)
      {
         i.test <- order(-y[(nTrain+1):cRows]) + nTrain
      }
      else
      {
         i.test <- NULL
      }
      i.timeorder <- c(i.train,i.test)

      y <- y[i.timeorder]
      Misc <- Misc[i.timeorder]
      x <- x[i.timeorder,,drop=FALSE]
      w <- w[i.timeorder]
      if(!is.na(offset)) offset <- offset[i.timeorder]
   }

   # create index upfront... subtract one for 0 based order
   if(ncol(x) > 1)
   {
      x.order <- apply(x[1:nTrain,],2,order,na.last=FALSE)-1
   }
   else
   {
      x.order <- order(x[1:nTrain,],na.last=FALSE)-1
   }
   x <- as.vector(data.matrix(x))
   predF <- rep(0,length(y))
   train.error <- rep(0,n.trees)
   valid.error <- rep(0,n.trees)
   oobag.improve <- rep(0,n.trees)

   if(is.null(var.monotone)) var.monotone <- rep(0,cCols)
   else if(length(var.monotone)!=cCols)
   {
      stop("Length of var.monotone != number of predictors")
   }
   else if(!all(is.element(var.monotone,-1:1)))
   {
      stop("var.monotone must be -1, 0, or 1")
   }
   fError <- FALSE
   gbm.obj <- .Call("gbm",
                    Y=as.double(y),
                    Offset=as.double(offset),
                    X=as.double(x),
                    X.order=as.integer(x.order),
                    weights=as.double(w),
                    Misc=as.double(Misc),
                    cRows=as.integer(cRows),
                    cCols=as.integer(cCols),
                    var.type=as.integer(var.type),
                    var.monotone=as.integer(var.monotone),
                    distribution=as.character(distribution),
                    n.trees=as.integer(n.trees),
                    interaction.depth=as.integer(interaction.depth),
                    n.minobsinnode=as.integer(n.minobsinnode),
                    shrinkage=as.double(shrinkage),
                    bag.fraction=as.double(bag.fraction),
                    nTrain=as.integer(nTrain),
                    fit.old=as.double(NA),
                    n.cat.splits.old=as.integer(0),
                    n.trees.old=as.integer(0),
                    verbose=as.integer(verbose),
                    PACKAGE = "gbmplus")
   names(gbm.obj) <- c("initF","fit","train.error","valid.error",
                       "oobag.improve","trees","c.splits")

   gbm.obj$bag.fraction <- bag.fraction
   gbm.obj$distribution <- distribution
   gbm.obj$interaction.depth <- interaction.depth
   gbm.obj$n.minobsinnode <- n.minobsinnode
   gbm.obj$n.trees <- length(gbm.obj$trees)
   gbm.obj$nTrain <- nTrain
   gbm.obj$response.name <- response.name
   gbm.obj$shrinkage <- shrinkage
   gbm.obj$train.fraction <- train.fraction
   gbm.obj$var.levels <- var.levels
   gbm.obj$var.monotone <- var.monotone
   gbm.obj$var.names <- var.names
   gbm.obj$var.type <- var.type
   gbm.obj$verbose <- verbose
   gbm.obj$Terms <- NULL

   if(distribution == "coxph")
   {
      gbm.obj$fit[i.timeorder] <- gbm.obj$fit
   }

   if(keep.data)
   {
      if(distribution == "coxph")
      {
         # put the observations back in order
         gbm.obj$data <- list(y=y,x=x,x.order=x.order,offset=offset,Misc=Misc,w=w,
                              i.timeorder=i.timeorder)
      }
      else
      {
         gbm.obj$data <- list(y=y,x=x,x.order=x.order,offset=offset,Misc=Misc,w=w)
      }
   }
   else
   {
      gbm.obj$data <- NULL
   }

   class(gbm.obj) <- "gbm"
   return(gbm.obj)
}

