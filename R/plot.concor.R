#' @title Plot Concordances
#'
#' @description This function plot the \code{ref} tree with the percents or counts of
#' concordances as node labels.
#'
#' @param x an object from "\code{concor}".
#' @param levelcol a character vector specifying the colors for node labels showing different
#'     levels of concordances. Defaults as \code{NULL}.
#' @param thres a numeric vector specifying the levels for node labels to be plotted with
#'     different colors. Must be consistent with \code{type}. Defaults as \code{NULL}.
#' @param type a character specifying the type of node labels. Can be either "percent" or "count".
#' @param node.cex a numeric specifying the scale of node labels. Defaults to \code{0.5}.
#' @param node.adj a numeric vector specifying the position of node labels. Defaults
#'     to \code{c(0.5, 0.5)}. Also see \code{plot.phylo}.
#' @param node.frame a character specifying the type of frame of node labels. Must be one of
#'     "rect", "circle", "none", or unambiguous abbreviation. Also see \code{plot.phylo}.
#'     Defaults to "\code{rect}".
#' @param ... other arguments to be passed to \code{plot} or to \code{plot.phylo}.
#' @export
#' @method plot concor
#' @importFrom ape plot.phylo
#' @importFrom ape nodelabels
#' @importFrom RColorBrewer brewer.pal
#' @importFrom methods hasArg
#' @importFrom graphics title
#' @details \code{thres} is required with \code{levelcol} arg. On the other hand, \code{levelcol}
#'     can left as default when \code{thres} is set, then \code{brewer.pal} (in \code{RColorBrewer}
#'     pcakage) will be called to produce \code{levelcol}. Length of \code{thres} must be one
#'     more than \code{levelcol}.
#' @seealso \code{\link{concor.node}}, \code{\link[ape]{plot.phylo}}


# bg will not be functional when thres or levelcol is provided

plot.concor <- function(x, levelcol = NULL, thres = NULL, type = c("percent", "count"),
        node.cex = 0.5, node.adj = c(0.5, 0.5), node.frame = "rect", ...){

  # removing arguments: "digit = 2", "title = TRUE"
  # library(ape)
  # suggest(RColorBrewer)

  if(!inherits(x, "concor"))
    stop("x must be an object of class \"concor\".")

  tree <- x$refTree
  datas <- x$data

  type <- match.arg(type, choices = c("percent", "count"), several.ok = FALSE)

  levels <- NULL
  if(methods::hasArg(levelcol)){
    # if levelcol is provided, thres is required
    lname <- names(levelcol)

    if(is.null(lname)){
      levels <- c(1:length(levelcol))
    } else{levels <- lname}

    if(is.null(thres))
      stop("thres must be provided with levelcol.")
    if(!inherits(thres, "numeric"))
      stop("thres must be numeric.")
    if(length(levels)-1 != length(thres))
      stop("Length of levelcol must be one more than that of thres.")
  }

  if(methods::hasArg(thres)){
    if(is.null(levelcol)){
      if(length(thres) < 3){
        levelcol <- RColorBrewer::brewer.pal(3, "GnBu")
        nn <- 3-length(thres)+1;nn
        nn <- c(nn:3);nn
        levelcol <- levelcol[nn]
      } else {
        levelcol <- RColorBrewer::brewer.pal(length(thres), "GnBu")
      }
      levelcol <- c("white", levelcol)
    }

    oo <- NULL
    if(length(thres) > 1){
      oo <- order(thres, decreasing = TRUE)
      thres <- thres[oo]
    }

    if(oo[1] > oo[2]){
      levelcol <- levelcol[c(length(levelcol):1)]
    }
  }

  # to set args for nodelabels
  if(methods::hasArg(bg)) bg <- list(...)$bg
  else bg <- "yellow"
  node.frame <- match.arg(node.frame, c("rect", "circle", "none"))
#  if(node.frame == "circle")
#    node.pch <- 21

  titles <- names(datas)
  for(d in 1:length(datas)){
    data <- datas[[d]]

    if(type == "percent"){
      per <- round(data$percent, digits = 2)
      data$n.concor <- per
      data <- data[ ,c(1:2)]
    }
    if(type == "count"){
      count <- data[ ,2]
      deno <- data[ ,3]
      lab <- paste(count, deno, sep = "/")
      data <- data.frame(data[ ,c(1:2)], lab)
    }
#    data <- data[ ,c(1:2)]

    if(!is.null(levelcol)){
      listl <- length(levelcol)
    } else {
      listl <- 1
    }

    ntext.list <- vector("list", listl)
    # from high score to low
    pthres <- NULL
    for(t in 1:length(thres)){
      cthres <- thres[t]
      if(is.null(pthres)){
        nt <- data[data[ ,2] >= cthres, ]
      } else {
        nt <- data[data[ ,2] >= cthres & data[ ,2] < pthres, ]
      } 
      
      if(inherits(nt, "numeric")){
        nt <- t(as.matrix(nt))
      }
      if(nrow(nt) != 0){
        ntext.list[[t]] <- nt
      }

      pthres <- cthres
    }
    # then for the final one levelcol (because thres is less than levelcol in length 1)
    cthres <- 0
    if(is.null(pthres)){
      nt <- data[data[ ,2] >= cthres, ]
    } else {
      nt <- data[data[ ,2] >= cthres & data[ ,2] < pthres, ]
    }

    if(inherits(nt, "numeric")){
      nt <- t(as.matrix(nt))
    }
    if(nrow(nt) != 0){
      ntext.list[[listl]] <- nt
    }

    # to plot
    ape::plot.phylo(tree, ...)
    if(hasArg(thres)){
      for(t in 1:length(ntext.list)){
        col <- levelcol[t]
        ntext <- ntext.list[[t]]
        if(!is.null(ntext)){
          if(type == "percent"){
            ape::nodelabels(text = ntext[ ,2],
                       node = ntext[ ,1], adj = node.adj,
                       frame = node.frame, bg = col, cex = node.cex)
          }
          if(type == "count"){
            ape::nodelabels(text = ntext[ ,3],
                       node = ntext[ ,1], adj = node.adj,
                       frame = node.frame, bg = col, cex = node.cex)
          }
        }
      }
    } else { # if no levelcol and thres, then plot all node in the same way
      ntext <- ntext.list[[1]]
      if(type == "percent"){
        ape::nodelabels(text = ntext[ ,2],
                   node = ntext[ ,1], adj = node.adj,
                   frame = node.frame, bg = bg, cex = node.cex)
      }
      if(type == "count"){
        ape::nodelabels(text = ntext[ ,3],
                   node = ntext[ ,1], adj = node.adj,
                   frame = node.frame, bg = bg, cex = node.cex)
      }
    }

#    if(title)
    graphics::title(main = paste(titles[d], " BP", sep = ""))
  }
}
