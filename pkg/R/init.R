.onAttach <- function(...) {
  welcome <- paste(""                                              ,
                   "----------------------------------------------",
                   "  'FusedANOVA' package version 0.00           ",
                   ""                                              ,
                   " Still under development... feedback welcome  ",
                   "----------------------------------------------",
                   sep = "\n")
  packageStartupMessage(welcome)
}

