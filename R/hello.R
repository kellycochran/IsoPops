.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste("Welcome to IsoDB version ",
                              utils::packageVersion("IsoDB"), ".", sep = ""))
}
