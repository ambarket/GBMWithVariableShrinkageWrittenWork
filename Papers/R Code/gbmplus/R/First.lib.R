".First.lib" <-
function(lib, pkg)
{
     library.dynam("gbmplus", pkg, lib)
     require(stats)
     require(lattice)
}

