pkg.pepit <- new.env()

pkg.pepit$RESIDUES=""
pkg.pepit$CONTACT=5.0
pkg.pepit$PRECISION=1.0
pkg.pepit$ADD="calpha"
pkg.pepit$ACC=10
pkg.pepit$POSE=TRUE
pkg.pepit$NBCLIQUES=5
pkg.pepit$NBHITS=20
pkg.pepit$SCORE=mapping_dist_sum2
pkg.pepit$MAXDELTADIST=25
pkg.pepit$MINDIST=0
pkg.pepit$MINCLIQUE=4
pkg.pepit$BCMIN=0
pkg.pepit$INTERCLIQUE=4
pkg.pepit$NBNEI=4
pkg.pepit$MAXEDGES=3000000
pkg.pepit$TYPES=c("a","b","c","o","n","A","C","O","N")
pkg.pepit$MINSCORE = 10
pkg.pepit$MAXCLASHES = 10
pkg.pepit$radius=1.4
pkg.pepit$PVALUE=FALSE
pkg.pepit$AAGAP=10
## utile ?
EMPTY_GRAPH=-3
NO_CLIQUE=-2
TOO_BIG=-1


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
