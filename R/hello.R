.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("Welcome to IsoPops version ",
                              utils::packageVersion("IsoPops"), ".", sep = ""))
}
