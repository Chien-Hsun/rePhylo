#' @title Summarize \code{cladeFilter} Results with Multiple Groupings
#'
#' @description This function summarizes the scenarios of all the groupings for each of
#' the \code{trees}, and generates barplots showing the numbers of removing sequences
#' by species (datasets), genetrees (genes) and groupings into pdf files.
#'
#' @param x an object from "\code{cladeFilter}".
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of test trees such as single gene trees.
#' @param plot.tree a logical specifying whether to plot \code{trees} with coloring tips
#'     to show tips suggested to be removed by \code{cladeFilter} for all groupings.
#'     Defaults to \code{TRUE}.
#' @param write.table a logical specifying whether to write out text files
#'     containing the summarized number and names of sequences suggested to be removed
#'     by \code{cladeFilter} for all groupings. Defaults to \code{TRUE}.
#' @param rm.col a character vector specifying the colors for tips removing by one grouping
#'     (defalts to \code{red}) or by multiple groupings (defaults to \code{blue}).
#' @param ... other arguments to be passed to \code{plot} or to \code{plot.phylo}.
#' @export
#' @importFrom ape plot.phylo
#' @importFrom ape tiplabels
#' @importFrom methods hasArg
#' @importFrom utils write.table
#' @importFrom graphics title
#' @seealso \code{\link{cladeFilter}}, \code{\link{barplotMulticf}}
#' @examples
#' \dontrun{
#' data(Brassidata)
#' trees <- Brassidata$trees
#' taxa <- Brassidata$taxaTable
#' # perform cladeFilter analysis
#' ures <- cladeFilter(trees = trees, taxa = taxa, level = 2)
#' # summarize cladeFilter results with multiple groupings
#' usum <- summarizeMulticf(ures, trees, plot.tree = FALSE, write.table = FALSE)
#' # to make barplots
#' utest <- barplotMulticf(x = usum, reorder="up", plot = TRUE, pdf=FALSE)
#' }




# for multiple constraints (groupings)
# there might be duplicated occurrences for some tips in different constraints
# either both removed or just removed by one constraint
# so the function here is to summarize from multiple constraints

summarizeMulticf <- function(x, trees, plot.tree = TRUE, write.table = TRUE,
            rm.col = c(rmOnce = "red", rmMulti = "blue"), ...){

  # library(ape)

  # x == res (class: list), contains:
  # lists: groups (class: cladeFilter)
  # lists: gene trees
  # lists: 1."result(untar_list)", 2."removingTips", 3."remainingTips",
  #         4."mrcaInTree", 5."node(fnode)", 6."subtree"
  # in "untar_list": a list named as targetNode
  # 1."num_wanted", 2."removeTips", 3."num_of_rm", 4."num_of_benefit"

  clas <- lapply(trees, function(z) inherits(z, "phylo"))
  clas <- unlist(clas)
  if(!inherits(trees, "list") | !all(clas))
    stop("Trees should be a list of objects of class \"phylo\".")

  filenames <- names(trees)
  if(is.null(filenames)){
#    stop("Trees should be named")
    # write in description
    names(trees) <- c(1:length(trees))
    filenames <- names(trees)
  }

  clas <- lapply(x, function(z) inherits(z, "cladeFilter"))
  clas <- unlist(clas)
  if(!inherits(x, "list") | !all(clas))
    stop("Input should be a list of objects of class \"cladeFilter\".")


  all.rm.tips <- lapply(x, function(z) .combine.(z, n = 2))
  check.rm <- unlist(lapply(all.rm.tips, function(z) lapply(z, is.null)))
  if(all(check.rm)){
    print("THERE IS NO TIP TO BE REMOVED.")
    return(NULL)
  } else {
    wanted_list <- lapply(x, function(z) .combine.(z, n = 3))
#    x == res??

    allrm<-NULL
    gnames <- names(all.rm.tips)
    for(i in 1:length(all.rm.tips)){
      t <- all.rm.tips[[i]]
      tnames <- names(t)
      gname <- gnames[i]
      for(j in 1:length(t)){
        tname <- tnames[j]
        tt <- t[[j]]
        if(!is.null(tt)){
          temp <- cbind(gname, tname, tt)
          if(is.null(allrm)){
            allrm<-temp
          } else {
            allrm <- rbind(allrm, temp)
          }
        }
      }
    }
#    if(!is.null(allrm)){
      # if "is.null(allrm) == TRUE",
      # there is no any tip to be removed in all the results
      if(inherits(allrm, "character"))
        allrm <- t(matrix(allrm))
      allrm <- data.frame(groupings = allrm[ ,1], genetrees = allrm[ ,2], species = allrm[ ,3],
                          stringsAsFactors = TRUE); nrow(allrm)
#    }
######################################################

    uniseq <- paste(allrm[ ,2], allrm[ ,3], sep = "_")
    wdup <- duplicated(uniseq)
    all.dedup <- allrm[!wdup, ]; nrow(all.dedup)

    # summarize numbers for genetrees
    gtrees <- all.dedup[!duplicated(all.dedup[ ,2]),2]
    gtrees <- as.list(gtrees)
    n.trees <- lapply(gtrees, function(zz) .sum.num(zz, n = 2, all.dedup))
    n.trees <- unlist(n.trees)
    gtrees <- unlist(gtrees)
    tab.trees <- data.frame(genetrees = gtrees, rm.num = n.trees, stringsAsFactors = FALSE)

    tab <- all.dedup[ ,c(2:3)]
    if(write.table){
      utils::write.table(tab, "summarize_rmTips.out", quote = FALSE,
                  row.names = FALSE, col.names = FALSE, sep = " ")
      utils::write.table(tab.trees, "summarize_rmNum.out", quote = FALSE,
                  row.names = FALSE, col.names = FALSE, sep = " ")
    }
  } # end of else for if(all(check.rm))


  if(plot.tree){
    # plot to see the summarized removed tips in one pdf file
    # plot subtree and wanted/remove tips
    # red are non-duplicated, and blue indicate duplicated removing tips

    # preparing options for nodelabels
    if(methods::hasArg(frame)) frame <- list(...)$frame
    else frame <- "none"
    if(methods::hasArg(bg)) bg <- list(...)$bg
    else bg <- "white"
    if(methods::hasArg(offset)) offset <- list(...)$offset
    else offset <- 0

    if(methods::hasArg(adj)) adj <- list(...)$adj
    else adj <- 0
    if(length(adj) == 1){
      treeadj <- adj
      tipadj <- c(0, 0.5)
    }
    if(length(adj) == 2){
      tipadj <- adj
      treeadj <- 0
    }

#    op <- list(...)
#    op$adj <- NULL
#    op$frame <- NULL
#    op$bg <- NULL
#    op$offset <- NULL

#    ... <- unlist(op)

    # preparing tip color
    if(is.null(names(rm.col))){
      rmOnce <- rm.col[1]
      rmMulti <- rm.col[2]
    } else {
      rmOnce <- rm.col[names(rm.col) == "rmOnce"]
      rmMulti <- rm.col[names(rm.col) == "rmMulti"]
    }

    dup <- allrm[wdup, ]
    dedup <- allrm[!wdup, ] # this will still contain the same ones in dup (only remove the redundant)
    duptip <- paste(dup[ ,2], dup[ ,3], sep = "")
    deduptip <- paste(dedup[ ,2], dedup[ ,3], sep = "")
    which <- which(deduptip %in% duptip) # change the object "duptip" to "which"
    if(length(which) > 0)
      dedup <- dedup[-which, ]

    for(i in 1:length(filenames)){
      # will plot all trees, with removing tips as red and blue
      # red: remove in one grouping; blue: remove in more than one groupings

      tree <- trees[[i]]
      tips <- tree$tip.label
      gt <- filenames[i]
      duptip <- dup[dup[ ,2] == gt,3]
      deduptip <- dedup[dedup[ ,2] == gt,3]

      col = rep("black", length(tips))

      which <- which(tips %in% deduptip)
      if(length(which) > 0){
        col[which] <- rmOnce
      }
      which <- which(tips %in% duptip)
      if(length(which) > 0){
        col[which] <- rmMulti
      }

      ape::plot.phylo(tree, tip.color = "white", adj = treeadj, ...)
      ape::tiplabels(col = col, text = tips, frame = frame, adj  = tipadj, bg = bg, offset = offset)
      graphics::title(main = gt)
    }
  } # end of plot.tree

  # to make summary barplots
  #2019.03.24 decide to make bar plots as an independent function
 # barplotMulticf(all.dedup, plot = TRUE, pdf = TRUE)

  return(all.dedup)

} # end of function




#####################################################################
# internal function
#' @title .combine.
#' @description Internal function of \code{youshu}
#' @param z z
#' @param n n
#' @export
.combine. <- function(z, n){
  t <- vector("list", length(z))
  names(t) <- names(z)
  for(i in 1:length(z)){
    zz <- z[[i]]; zz
    for(j in 1:length(zz)){
      zzz <- zz[[n]]
      if(!is.null(zzz)){
        t[[i]] <- zzz
      }
    }
  }
  return(t)
}
