.onLoad <- function(libname, pkgname){
  eaconcache.pkgname <- 'EaCoN.cache'
  if(eaconcache.pkgname %in% installed.packages()) {
    eaconcache.dir <- system.file(package = eaconcache.pkgname)
    cache.dir <- paste0(eaconcache.dir, '/extdata/')
    options(EACON_CACHE=cache.dir)
    packageStartupMessage("Welcome to EaCoNic !\n*-*-*-*\n\nThis package requires a cache of chromosome and genome structure, as well as genomic annotations to be used.\nThis is available for most of UCSC genomes thanks to the 'EaCoN.cache package.\nSee https://github.com/GustaveRoussy/EaCoN.cache/" )
  } else {
    
  }
}