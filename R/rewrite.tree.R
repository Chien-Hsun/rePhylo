#' @title rewrite.tree
#' @description Internal function of \code{youshu}. Make the numbering of nodes and tips of a reordered tree
#'     as re-read in.
#' @param tree tree
#' @export
#' @keywords internal
#' @seealso \code{\link{root_dist}}, \code{\link[ape]{root.phylo}}

# to get the post ordered tree without write out a file
# applicable for root, ladderize, rotateConstr, ...etc

rewrite.tree <- function(tree){
  after_tree <- tree

  is_tip <- tree$edge[ ,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip,2]
  tips2 <- tree$tip.label[ordered_tips]
  # now the tips2 is the new order of tips after rooting and writing/reading

  edge2 <- tree$edge

  res <- matrix(nrow = nrow(edge2), ncol = 2)
  # to change the order of internal nodes
  e1 <- edge2[ ,1]
  dd <- e1[!duplicated(e1)] ## old e1
  ne1 <- c(1:length(dd))
  ne1 <- ne1+dd[1]-1 ## new e1
  #  e2 <- edge2[ ,2] ## old e2
  res1 <- rep(NA, length(e1))

  turn <- function(x, res1, e1, dd, ne1){
    w <- which(e1 == dd[x])
    res1[w] <- ne1[x]
    return(res1)
  }
  n <- c(1:length(dd))
  temp <- lapply(as.list(n), function(x) turn(x, res1, e1, dd, ne1))

  newe <- c(0, 0)
  for(i in 1:length(temp)){
    tt <- temp[[i]]
    w <- which(!is.na(tt))
    ttt <- tt[w]
    dat <- data.frame(w, ttt)
    newe <- rbind(newe, dat)
  }
  newe <- newe[-1, ]
  ordere <- order(newe[ ,1])
  newee <- newe[ordere, ]
  # now the newee[ ,2](ttt) == edge3[ ,1]
  res[ ,1] <- newee[ ,2]

  # for internal nodes in e2
  e2 <- edge2[ ,2]
  for(i in 1:nrow(newee)){
    old <- e1[newee[i,1]]
    new <- newee[i,2]
    w <- which(e2 == old)
    res[w,2] <- new
  }

  # for tips
  e2 <- edge2[ ,2]
  et2 <- e2[e2 <= length(tips2)] # old e2
  ne2 <- c(1:length(et2)) ## new e2
  n <- c(1:length(et2))
  res1 <- rep(NA, length(e1))
  temp <- lapply(as.list(n), function(x) turn(x, res1, e2, et2, ne2))
  newe <- c(0, 0)
  for(i in 1:length(temp)){
    tt <- temp[[i]]
    w <- which(!is.na(tt))
    ttt <- tt[w]
    dat <- data.frame(w, ttt)
    newe <- rbind(newe, dat)
  }
  newe <- newe[-1, ]
  ordere <- order(newe[ ,1])
  newee <- newe[ordere, ]
  for(i in 1:nrow(newee)){
    w <- newee[i,1]
    new <- newee[i,2]
    res[w,2] <- new
  }


  after_tree$edge <- res
  after_tree$tip.label <- tips2
  if("node.label" %in% names(after_tree)){
    # for node labels
    nl <- tree$node.label[c(2:length(tree$node.label))]
    e2 <- tree$edge[tree$edge[ ,2] > length(tree$tip.label),2]
    old <- order(e2)

    newnl <- rep(NA, length(nl))
    for(i in 1:length(nl)){
      ee <- old[i]
      newnl[ee] <- nl[i]
    }

    newnl <- c(tree$node.label[1], newnl)
    after_tree$node.label <- newnl

  }

  return(after_tree)
}


