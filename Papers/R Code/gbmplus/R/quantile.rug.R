"quantile.rug" <-
function(x,prob=(0:10)/10,...)
{
     quants <- quantile(x[!is.na(x)],prob=prob)
     if(length(unique(quants)) < length(prob))
     {
          quants <- jitter(quants)
     }
     rug(quants,...)
}

