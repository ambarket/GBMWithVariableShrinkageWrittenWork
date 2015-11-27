"gbm.loss" <-
function(y,f,w,offset,dist,baseline)
{
   if(!is.na(offset))
   {
      f <- offset+f
   }
   switch(dist,
          gaussian = weighted.mean((y - f)^2,w) - baseline,
          bernoulli = -2*weighted.mean(y*f - log(1+exp(f)),w) - baseline,
          laplace = weighted.mean(abs(y-f),w) - baseline,
          adaboost = weighted.mean(exp(-(2*y-1)*f),w) - baseline,
          poisson = -2*weighted.mean(y*f-exp(f)) - baseline,
          stop(paste("Distribution",dist,"is not yet supported for method=permutation.test.gbm")))
}

