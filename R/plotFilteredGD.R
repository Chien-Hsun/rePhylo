#' @title plot the gd result corrected by filter_GD_by_coverage
#'     and denominator
#'
#' @description \code{plotFilteredGD} plot a tree with node labels as
#'     the corrected values (both count and percent). Only nodes both 
#'     shown in both \code{gd} and \code{deno} have node labels. Only 
#'     nodes with deno > 0 are shown.
#'
#' @param tree the species tree used to do tree reconciliation. 
#' @param gd the summary table of \code{filter_GD_by_coverage}.
#' @param deno the result of \code{denominator}.
#' @param type can be one of "percent", "count", or "both", determines 
#'     the type of threshod value used for node color. Defaults to \code{"percent"}.
#' @param thres_node optional. A numeric specify which node 
#'     (in ape id) are used as threshold to give colors on node labels.
#' @param thres_count nodes with counts above this threshold will be plotted 
#'     as red. Defaults to \code{100}
#' @param thres_percent nodes with percent above this threshold will be plotted 
#'     as red. Defaults to \code{5}
#' @param color two characters specify the colors for node labels. The first is 
#'     for node above the threshold, and the second is for those below.
#' @param tip.offset offset of tip labels in \code{ggtree}. Defaults to \code{NULL}.
#' @param plot.margin four numeric indicate plot margin. Defaults to \code{c(6, 140, 6, 6)}. 
#' @param xlim used as the x scale in \code{ggtree}.
#' @param node.cex size of node labels.
#' @export
#' @importFrom ggtree ggtree
#' @importFrom ggtree geom_tiplab
#' @importFrom ggtree geom_label2
#' @importFrom ggtree theme_tree2
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 coord_cartesian


plotFilteredGD <- 
  function(tree, gd, deno,
           type = "percent",
           thres_node, thres_count = 100,
           thres_percent = 5,
           color = c("red", "yellow"),
           tip.offset = NULL,
           node.cex = 1,
           plot.margin = c(6, 140, 6, 6),
           xlim = NULL){
  
  cat("Making tree plot using ggtree.....\n")
  
  edge <- tree$edge
  edge <- edge[edge[ ,2] > length(tree$tip.label), 2]
  
  lab <- count <- percent <- rep(NA,length(edge))
  for(i in 1:length(edge)){
    n <- edge[i]
    c <- gd[which(gd[ ,"node"] == n),3]
    de <- deno[which(deno[ ,"node"] == n),2]
    if(length(c) > 0 & length(de) > 0){
      if(de > 0){
        p <- (c / de) * 100
        percent[i] <- p
        count[i] <- c
        lab[i]<-paste0(c,"/",de," [",round(p,3),"%]")
      }
    }
  }
  # add root node
  edge <- c(length(tree$tip.label)+1, edge)
  count <- c(NA, count)
  percent <-c(NA, percent)
  lab <- c(NA, lab)
  # a complete table
  tab <- data.frame(edge = edge,
                    count = count,
                    percent = percent,
                    lab = lab)

  # to plot using ggtree
  ll<-length(tree$tip.label)
  if(hasArg(thres_node)){
    thresc<-tab[which(tab[,1] == thres_node),"count"]
    thresp<-tab[which(tab[,1] == thres_node),"percent"]
  } else {
    thresc<-thres_count
    thresp<-thres_percent
  }
  if(type == "both"){
    whichred<-which(as.numeric(tab[,"count"]) >= as.numeric(thresc) & as.numeric(tab[,"percent"]) >= as.numeric(thresp))
  }
  if(type == "percent"){
    whichred<-which(as.numeric(tab[,"percent"]) > as.numeric(thresp))
  }
  if(type == "count"){
    whichred<-which(as.numeric(tab[,"count"]) >= as.numeric(thresc))
  }
  
  tree2<-tree
  
  whichred <- whichred + ll; whichred
  whichyel <- c(1:(tree2$Nnode + ll))
  whichyel <- whichyel[-whichred]; whichyel
  labyel <- labred <- rep(FALSE, (tree2$Nnode + ll))
  labred[whichred] <- TRUE; labred
  labyel[whichyel] <- TRUE; labyel
  labyel[c(1:ll)] <- FALSE # don't label tips

  tree2$node.label <- as.character(tab[ ,"lab"])
  
  if(is.null(tip.offset)){
    tip.offset <- 0
  }
  
  g <- ggtree::ggtree(tree2)
  gg <- g +
    ggtree::geom_tiplab(offset = tip.offset) +
    ggtree::geom_label2(ggplot2::aes(label = label, subset = labred), #label %in% label[grep(label,pattern="%")] & tfc == "TRUE" & tfp == "TRUE"),
                bg = color[1], size = node.cex) +
    ggtree::geom_label2(ggplot2::aes(label = label, subset = labyel), #label %in% label[grep(label,pattern="%")] & -whichc & -whichp),
                bg = color[2], size = node.cex) +
    ggplot2::coord_cartesian(clip = "off") + 
    ggtree::theme_tree2(plot.margin = ggtree::margin(plot.margin)) +
    ggplot2::labs(title = "corrected_w_sub_coverage") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  if(!is.null(xlim)){
    gg <- gg + ggplot2::xlim(xlim[1], xlim[2])
  }
  
  print(gg)

  return(tab)
}

