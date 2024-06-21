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
read.config <- function(f="pepit_cfg.R") {
  if (file.exists(f)) {
      source(f, local=pkg.pepit)
  }
}

pkg.pepit <- new.env()
read.config("pepit_cfg.R")
