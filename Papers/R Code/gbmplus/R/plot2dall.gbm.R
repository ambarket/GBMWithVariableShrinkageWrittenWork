"plot2dall.gbm" <-
function (obj, cols, box.add = TRUE, pts.add = TRUE, pty, ...) 
{
    if (!is.null(names(obj))) {
        nms <- obj$var.names
        fact <- obj$var.type > 0
    }
    else {
        nms <- obj[[1]]$var.names
        fact <- obj[[1]]$var.type > 0
    }
    if (missing(cols)) 
        cols <- 1:length(nms)
    if (missing(pty)) 
        pty <- c("c", "l")[fact+1]
    nc <- length(cols)
    if (nc == 2) 
        plot.gbm(obj, cols = cols, ...)
    else {
        old.par <- par(mfcol = c(nc - 1, nc - 1), oma = c(3, 
            3, 1, 1), mar = rep(0.5, 4))
        on.exit(par(old.par))
        for (i in 1:(nc - 1)) {
            if (i > 1) for (k in 1:(i-1)) frame()
            for (j in (i + 1):nc) {
                ncs <- c(i,j)
                ptij <- ifelse(any(pty[c(i, j)] == "l"), "l", 
                  pty[i])
                plot.gbm(obj, ncs, xlab = "", ylab = "", axes = F, 
                  pts = pts.add, pty = ptij, ...)
                if (box.add) 
                  box()
                if (i == 1 & !fact[j]) 
                  mtext(nms[j], 2, line = 1, cex = 1.25 * par()$cex)
                if (j == nc  & !fact[i]) 
                  mtext(nms[i], 1, line = 1, cex = 1.25 * par()$cex, 
                    las = 0)
            }
        }
      }
}

