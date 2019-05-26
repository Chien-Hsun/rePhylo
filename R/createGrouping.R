#' @title Creating a Taxa Table as Input to Other Functions in package \code{youshu}
#'
#' @description By giving \code{percent} or \code{count}, the groupings are created from nodes
#'     showing the percent or count of concordances above \code{percent} or \code{count},
#'     respectively, and returns a \code{data.frame} that can be used directly as \code{taxa}
#'     in functions \code{cladeFilter} and \code{backboneBP}. Also can specify the nodes to
#'     produce groupings by giving \code{nodes} directly.
#'
#' @param concor an object from "\code{concor.node}".
#' @param percent a numeric between 0 to 1. The groupings are created from nodes
#'     showing the percent of concordances above \code{percent}. Defaults to \code{0.5}.
#' @param count a numeric or integer higher than or equals to 1.
#'     The groupings are created from nodes showing the count of concordances
#'     above \code{count}. Defaults to \code{NULL} and uses \code{percent}.
#' @param bp a numeric specifying with which bp the result in \code{concor} is used
#'     to create \code{taxa}. Defaults to \code{NULL} and the first result in
#'     \code{concor} is used.
#' @param nodes a numeric vector specifying the node(s) to create \code{taxa},
#'     regardless of \code{percent} and \code{count}. Defaults to \code{NULL}.
#' @export
#' @importFrom ape extract.clade
#' @importFrom methods hasArg
#' @seealso \code{\link{concor.node}}
#' @examples
#' \dontrun{
#' data(Brassidata)
#' trees <- Brassidata$trees
#' ref <- Brassidata$ref
#' concor.n<-concor.node(ref = ref,trees = trees, bp = c(0,30,50), getTreeNames = FALSE)
#' # to get groupings from nodes with >70% tree in concordance
#' taxa <- createGrouping(concor.n, percent = 0.7)
#' }


createGrouping <- function(concor, percent = 0.5, count = NULL, bp = NULL, nodes = NULL){

  if(!inherits(concor, "concor"))
    stop("Input must be an object of class\"concor\".")

  data <- concor$data
  bps <- names(data)
  if(!methods::hasArg(bp)) data <- data[[1]]
  if(methods::hasArg(bp)){
    if(length(bp) > 1) bp <- bp[1]
    # write in description: if(length(bp) >1) bp <- bp[1]; if(!hasArg(bp)) data <- data[[1]]
    w <- which(bps == bp)
    if(length(w) == 0){
      cat("bp is not in the concor object. ", bps[1], " bp will be used.", sep = "")
      data <- data[[1]]
    }
    else {
      data <- data[[w]]
    }
  }
  tdata <- data

  if(!methods::hasArg(nodes)){
    if(!methods::hasArg(count) & !methods::hasArg(percent)){
      cat("No count or percent is provided. Percent = 0.5 will be used.")
      data <- data[data[ ,4] >= percent, ]
      data <- data[ ,c(1, 4)]
    }
    if(!methods::hasArg(count) && methods::hasArg(percent)){
      percent <- percent[1]
      data <- data[data[ ,4] >= percent, ]
      data <- data[ ,c(1, 4)]
    }
    if(methods::hasArg(count) && !methods::hasArg(percent)){
      count <- count[1]
      data <- data[data[ ,2] >= count, ]
      data <- data[ ,c(1, 2)]
    }
    if(methods::hasArg(count) && methods::hasArg(percent)){
      cat("Both count and percent is provided. Percent will be used.")
      percent <- percent[1]
      data <- data[data[ ,4] >= percent, ]
      data <- data[ ,c(1, 4)]
    }
  }

  # if hasArg(nodes)overwites the "percent" and "count" args and just choose these nodes
  if(methods::hasArg(nodes)){
    data <- tdata
    data <- data[data[ ,1] %in% nodes, ]
  }

  # to get the tips of each node
  tree <- concor$refTree
  tnodes <- as.list(data[ ,1])
  names(tnodes) <- tnodes
  getGroups <- function(x, tree){
    stree <- ape::extract.clade(phy = tree, node = x)
    tips <- stree$tip.label
    return(tips)
  }
  groups <- lapply(tnodes, function(x) getGroups(x, tree))
  # to combine all groups into one table
  # for preparation of taxa table
  tm <- c(0, 0)
  for(i in 1:length(groups)){
    tips <- groups[[i]]
    nn <- tnodes[[i]]
    gg <- rep(x = nn, length(tips))
    dat <- data.frame(tips = tips, group = gg, stringsAsFactors = FALSE)
    tm <- rbind(tm, dat)
  }
  tm <- tm[-1, ]
  tm
  return(tm)
} # end of function
