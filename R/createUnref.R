# currently not used in package (2018.08.05)

#' @title createUnref
#' @description Internal function of \code{BetterTree}.
#' @param concor an object of class "\code{concor.node}".
#' @param bp a numeric specifying the bp value of node(s) to be set as unresolved node
#'     when the corresponding node labels are lower than \code{bp}.
#' @param percent a numeric between 0 to 1. The constraints are created from nodes
#'     showing the percent of concordances above \code{percent}. Defaults to \code{0.5}.
#' @export
#' @importFrom methods hasArg
#' @importFrom ape di2multi

# create an unresolved tree as ref from the result of concor.node, to be used in cladeFilter
createUnref <- function(concor, bp, percent){
  if(!inherits(concor, "concor"))
    stop("Input must be an object of class\"concor\".")
  if(!methods::hasArg(percent))
    stop("percent is required.")
  if(percent > 1 | percent < 0)
    stop("percent must be a value within 0 and 1.")

  data <- concor$data
  if(hasArg(bp)){
    w <- which(names(data) == bp)
    if(length(w) > 0) data <- data[[w]]
    else data <- data[[1]]
  } else {
    data <- data[[1]]
  }
  per <- data$percent
  w <- which(per < percent)
  if(length(w) == 0)
    stop(paste("There is no node with percent higher than or equal to ", percent, ".", sep = ""))
  if(length(w) > 0){
    data <- data[w, ]
  }

  tree <- concor$refTree
  nodes <- data[ ,1]

  # to get edge lengths
  edgelen <- tree$edge.length
  edge <- tree$edge
  w <- which(edge[ ,2] %in% nodes)
  if(length(w) > 0)
    edgelen[w] <- 0
  tree$edge.length <- edgelen

  untree <- ape::di2multi(tree)
  return(untree)
}
