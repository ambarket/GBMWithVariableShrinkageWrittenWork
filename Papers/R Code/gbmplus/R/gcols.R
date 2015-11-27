"gcols" <-
function (x, FUNCOL = heat.colors(50)) 
{
    x <- x[-1, -1]
    dim(x) <- NULL
    x <- ceiling(length(FUNCOL) * (x - min(x, na.rm = T))/(max(x, na.rm = T) - 
        min(x, na.rm = T)))
    x[x == 0] <- 1
    FUNCOL[x]
}

