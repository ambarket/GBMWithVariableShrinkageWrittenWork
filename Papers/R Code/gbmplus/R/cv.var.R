"cv.var" <-
function(x, cv = 10) {
        x <-  match(x, sort(unique(x)))
        maxx <- max(x)
        if (maxx >= cv) {
            x <- match(x, sample(maxx))
            grps <- ceiling((cv * cumsum(table(x)))/length(x))
            x <- number(grps[x])
        }
        x
}

