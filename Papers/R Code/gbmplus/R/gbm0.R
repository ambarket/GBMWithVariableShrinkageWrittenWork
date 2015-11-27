"gbm0" <-
function(..., keep.data=TRUE, frac = 0.95, monitor=T) 
{
    obj <- gbm(..., keep.data = keep.data, monitor = FALSE)
    done <- (obj[[1]]$best/obj[[1]]$n.tree < frac)
    i <- 1
    cat("# Loops =",i,"")
    while (!done) {
    i <- i+1
    cat(i,"")
    obj <- gbm.more(obj)
    done <- (obj[[1]]$best/obj[[1]]$n.tree < frac) 
    }
    if (monitor){
        cat("\nN trees =", obj[[1]]$n.trees, ": Depth =", obj[[1]]$interaction.depth, 
            ": Minimum node size =", obj[[1]]$n.minobsinnode, "\n")
        cat("Distribution =", obj[[1]]$distribution, ": Shrinkage =", 
            signif(obj[[1]]$shrinkage), ": Bag Fraction =", obj[[1]]$bag.fraction, 
            "\n")
        cat("Best size =", obj[[1]]$best, ": Best error =", signif(obj[[1]]$best.error, 
            5), "\n\n")
        monitor.gbm(obj)
    }
    obj
}

