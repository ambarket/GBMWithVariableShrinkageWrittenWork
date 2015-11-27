"nbig.plt" <-
function (n) 
{
    if (n == 1) 
        m <- c(1, 1)
    else if (n == 2) 
        m <- c(1, 2)
    else if (n <= 4) 
        m <- c(2, 2)
    else if (n <= 6) 
        m <- c(2, 3)
    else if (n <= 9) 
        m <- c(3, 3)
    else m <- c(3, 4)
    m
}

