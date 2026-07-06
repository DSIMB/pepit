#' set a pepit parameter accessible for users
#'
#' @param tag a string : parameter name
#' @param value values set for the tag parameter
#' @description 
#' @return a value
#' @export
#'
#' @examples
set.pepit <- function(tag, value) {
  assign(tag , value, envir=pkg.pepit)
}

#' get a pepit parameter
#'
#' @param tag a string : parameter name
#' @description 
#' @return a value
#' @export
#'
#' @examples
get.pepit <- function(tag) {
  get(tag, envir=pkg.pepit)
}

#' read a config file 
#'
#' @param a config R file
#' @description side effect update environment pkg.pepit
#' @return nothing
#' @export
#'
#' @examples
read.config <- function(f) {
  if (file.exists(f)) {
      message("load config:",f,"\n")
      source(f, local=pkg.pepit)
  } else {
    message("No config file", f, "! default config loaded\n")
  }
}

pkg.pepit <- new.env()
