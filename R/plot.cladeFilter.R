#' @title Plot \code{cladeFilter} Result in \code{trees} for One Grouping
#'
#' @description Plot \code{trees} with tips colored for members of a grouping suggested to be
#' remained (green) or to be removed (red).
#'
#' @param x an object from "\code{cladeFilter}".
#' @param node.adj a numeric vector specifying the \code{adj} of node labels. Defaults
#'     to \code{c(0.5, 0.5)}. Also see \code{plot.phylo}.
#' @param ... other arguments to be passed to \code{plot} or to \code{plot.phylo}.
#' @export
#' @method plot cladeFilter
#' @importFrom ape plot.phylo
#' @importFrom methods hasArg
#' @importFrom graphics title
#' @seealso \code{\link{cladeFilter}}, \code{\link[ape]{plot.phylo}}

# to plot trees after cladeFilter process
# input "x" as an object of "cladeFilter" but not "a list of cladeFilter"
# ... will be passed to plot.phylo
plot.cladeFilter <- function(x, node.adj = c(0.5, 0.5), ...){

  # removing arguments: "tip.col = c(remain = "green", rm = "red")",
  # node.label = c("subclade", "all", "none")", "titles = NULL"

  # input as res, and the x == res

  if(!inherits(x, "cladeFilter"))
    stop("x should be a list of objects of class \"cladeFilter\".")

#  if(is.null(titles))
   titles <- names(x)

#  if(!hasArg(node.label))
#    node.label <- "subclade"
#  node.label <- match.arg(node.label, choices = c("subclade", "all", "none"))

  # preparing options for nodelabels
  if(methods::hasArg(frame)) frame <- list(...)$frame
  else frame <- "c"
  if(methods::hasArg(bg)) bg <- list(...)$bg
  else bg <- "yellow"


  # preparing tip color
#  if(is.null(names(tip.col))){
#    remain <- tip.col[1]
#    rm <- tip.col[2]
#  } else {
  tip.col = c(remain = "green", rm = "red")
  remain <- tip.col["remain"]
  rm <- tip.col["rm"]
#  }

  for(set in 1:length(x)){
    z <- x[[set]]
    untar_list <- z[[1]]
    removeTips <- z[[2]]
    wantedTips <- z[[3]]
    mrcaInTree <- z[[4]]
    node <- z[[5]]
    subtree <- z[[6]]
    if(!is.null(untar_list)){
      targetNode <- names(untar_list)
      targetNode <- as.numeric(targetNode)
    }

    if(!is.null(subtree)){
      # if subtree == NULL indicate zero or only 1 tip in tree, and no need to plot
      if(is.null(untar_list)){
        # if no untar_list but have subtree, indicates there are target tips  but nothing to be removed
        # then just plot subtree
        which <- which(subtree$tip.label %in% mrcaInTree)
        col <- rep("black", length(subtree$tip.label))
        col[which] <- remain
        ape::plot.phylo(subtree, tip.color  = col, ...)
#        if(node.label == "all")
#          ape::nodelabels(frame = frame, bg = bg, adj = node.adj)
#        if(node.label == "subclade")
#          if(!is.null(untar_list))
#            ape::nodelabels(node = targetNode, frame = frame, bg = bg, adj = node.adj)
        graphics::title(main = titles[set])

      } else {
        green <- which(subtree$tip.label %in% wantedTips)
        red <- which(subtree$tip.label %in% removeTips)
        mm <- match(green, red)
        if(!all(is.na(mm)))  # This should not happen
          warning(paste(set, ": Conflict tip occurs in both removedTips and wantedTips", sep = ""))
        col <- rep("black", length(subtree$tip.label))
        if(length(green) > 0) col[green] <- remain
        if(length(red) > 0) col[red] <- rm

        # plot subtree and wanted/remove tips
        ape::plot.phylo(subtree, tip.color = col, ...)
#        if(node.label == "all")
#          ape::nodelabels(frame = frame, bg = bg, adj = node.adj)
#        if(node.label == "subclade")
#          ape::nodelabels(node = targetNode, frame = frame, bg = bg, adj = node.adj)
        graphics::title(main = titles[set])

      }
    }
  }
}
