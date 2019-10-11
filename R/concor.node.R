#' @title Concordances of Trees with All Nodes of A Reference Tree
#'
#' @description \code{concor.node} analyzes and returns for each node on the ref tree
#' the number and names of \code{trees} showing concordances.
#'
#' @param ref an object of class "\code{phylo}". A reference tree.
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of test trees such as single gene trees.
#' @param bp a bootstrape value, which the concordance is counted
#'     when the corresponding node label of \code{trees} is higher than.
#' @param getTreeNames a logical specifying whether to return the names of \code{trees}
#'     showing concordances for each node on \code{ref}. Defaults to \code{FALSE}.
#' @param node a numeric or integer specifying the node(s) on \code{ref}
#'     to analyze the concordance. Defaults to \code{NULL} and all the nodes on
#'     \code{ref} are analyzed.
#' @export
#' @details Concordance for a node is defined as all descendent tips of the node of \code{ref}
#'     are grouped together without other tips in \code{trees}. The main analysis of this
#'     function requires rooted trees; when the \code{trees} are not rooted, \code{concor.node}
#'     roots \code{trees} according to \code{ref} during the analysis. \code{concor.node}
#'     returns an object of class "\code{concor}" as a list of: (1) \code{tree}: the \code{ref};
#'     (2) data: a \code{data.frame} showing the counts and percents of concordances for each
#'     node; (3) \code{conTrees}: names of \code{trees} that are concordant with \code{ref} for
#'     each node (when \code{getTreeNames = TRUE}).
#' @importFrom methods hasArg
#' @examples
#' data(Brassidata)
#' trees <- Brassidata$trees
#' ref <- Brassidata$ref
#' concor.n<-concor.node(ref = ref,trees = trees, bp = c(0,30,50), getTreeNames = FALSE)



# 0803: NEED to check with: trees with no bp values

###################################################################################
# This part is just set the "default" function of concor.node to be as two steps of rooting version

concor.node <- function(ref, trees, bp = 0, getTreeNames = FALSE, node = NULL){
# removing "warn = FALSE", "unresolve.ref = FALSE", "unresolve.para = c(NA,0.8)"

  clas <- lapply(trees, function(z) inherits(z, "phylo"))
  clas <- unlist(clas)
  if(!all(clas))
    stop("Trees should be a list of objects of class \"phylo\".")

  if(!inherits(ref, "phylo"))
    stop("Ref tree should be an objects of class \"phylo\".")

  treenames <- names(trees)
  if(is.null(treenames)){
    treenames <- c(1:length(trees))
  }

  rdata <- .concor.node.unrooted(ref = ref, trees = trees, bp = bp, getTreeNames = getTreeNames,
                                node = node, warn = TRUE,
                                unresolve.ref = FALSE, unresolve.para = NULL,
                                from = "concor.node")

  class(rdata) <- "concor"
  return(rdata)
} # end of function



#####################################################################################
# one ref to multiple trees
# if you want multi v.s. multi, just do lapply() or for() for ref trees
# "getTreeNames" can return the object with which trees are in concordance with ref tree for each node
# "node" can specify which node to calculate the concordance (if NULL, check the whole tree)

# .concor.node.unrooted: (conparing to concor.node.rooted)
# add two steps of rooting, with .denominator function
# and, to judge concordances based on the trees in dres (.denominator pools)
# and, change the selection of new roots from the descendent nodes of root of the ref tree
# (and decide the nodes to grep from trees2 based on the "grouping" of the two des nodeds)
#' @title .concor.node.unrooted
#' @description Internal function of \code{youshu}
#' @param ref ref
#' @param trees trees
#' @param bp bp
#' @param getTreeNames getTreeNames
#' @param node node
#' @param warn warn
#' @param unresolve.ref unresolve.ref
#' @param unresolve.para unresolve.para
#' @param from from
#' @importFrom methods hasArg
#' @export

.concor.node.unrooted <- function(ref, trees, bp = 0, getTreeNames = FALSE, node = NULL, warn = TRUE,
                               unresolve.ref = FALSE, unresolve.para = c(bp = NA, percent = 0.8),
                               from){
  # unresolve.para as c(bp, percent)

  ##############################################################################
  # check inputs

  clas <- lapply(trees, function(z) inherits(z, "phylo"))
  clas <- unlist(clas)
  if(!all(clas))
    stop("Trees should be a list of objects of class \"phylo\".")

  if(!inherits(ref, "phylo"))
    stop("Ref tree should be an objects of class \"phylo\".")

  treenames <- names(trees)
  if(is.null(treenames)){
    treenames <- c(1:length(trees))
  }

  # check if trees have bp values
  cbps <- lapply(trees, function(x) x$node.label)
  cbps <- unlist(lapply(cbps, length))
  if(any(cbps  ==  0)){
    cat("Note: Your test trees have no bp values. Arg \"bp\" is ignored.\n")
    #    ref$node.label <- rep(0, ref$Nnode)
    tt <- function(x){
      x$node.label <- rep(0, x$Nnode)
      return(x)
    }
    trees <- lapply(trees, function(x) tt(x))
    bp <- 0
  }
#  ref <- addbp(ref)

  # check is.rooted for ref
#  if(!is.rooted(ref)){
#    warning("ref tree is not rooted.")
#  }

  if(unresolve.ref){
    if(!methods::hasArg(unresolve.para)){
      cat("No arg \"unresolve.para\" is provided with unresolve.ref = TRUE. BP = ", bp[1], " and percent = 0.8 will be used.\n", sep = "")
      unbp <- bp[1]
    }
    if(methods::hasArg(unresolve.para)){
      unbp <- unresolve.para[1]
      if(!unbp %in% bp){
        cat("BP in unresolve.para is not one of the ", bp, ". BP = ", bp[1], " will be used.\n", sep = "")
        unbp <- bp[1]
      }
      percent <- unresolve.para[2]
    }
  }

  #################################################################################
  # check if tip labels are all matched
  temp <- .match.tiplabels(trees = trees, ref = ref)
  trees <- temp$trees
  ref <- temp$ref

  #################################################################################
  # main content

  refinfo <- .getInfo(ref, node = node)
  # info contains a list of two: matrix(nodeID+bp) and list (node/tipnames)

  tiplist <- refinfo[[2]]
  # tiplist is a list in length and IDs of all nodes of the original ref tree
  # and then need to remove tips lacking in the test trees

  setBPs <- bp
  reftips <- ref$tip.label

  conTreesl <- con_results <- data <- vector("list", length(setBPs))
  names(conTreesl) <- names(data) <- names(con_results) <- setBPs

  # to get the .denominator of ref nodes,
  # i.e. how many trees have tips that can define this node
  # rationale: below the target node, both of the the two child clades (tips) have at least one tip for each subclade...
  # regardless of structure, so even unrooted tree is fine
  dres <- lapply(trees, function(x) .denominator(x, ref, refinfo))
  # first level as trees, second as nodes

  # to get the number of denominators (the number of trees that have the "node" on ref)
  ttt <- function(x, dres){
    tres <- lapply(dres, function(xx) return(xx[x])) # xx[x] as trees[nodes]
    # tres as T/F of all trees for the same nodes
    tres <- unlist(tres)
    res <- length(which(tres))
    return(res)
  }
  len <- as.list(c(1:length(dres[[1]]))) #nodes
  ddres <- unlist(lapply(len, function(x) ttt(x, dres)))
  # ddres is now the denominators of nodes

  # to make a list contains tiplist and which trees are included as denominators
  allnodes <- names(tiplist)
  tiplist2 <- tiplist
  for(ll in 1:length(tiplist)){
#    nn <- allnodes[ll]

    temp <- unlist(lapply(dres, function(x) x[[ll]]))
    ww <- which(temp)
    # technically, length(tr) should not  ==  0, because we have perform match.tiplabels
    # so that for each tip in ref tree, there is at least one trees contain one of them.
    # But it is still possible for derived nodes (e.g. having only 2 tips) to have no trees having the representatives for same node
    # if the two tips lost in different trees

    if(length(ww)  ==  0){
      stop("None of the trees has representatives for node ", allnodes[ll], ".")
    } else {
      tiplist2[[ll]] <- list(tips = tiplist[[ll]], trees.w.node = ww)
    }
  }


  #################################################################################

    # for unrooted test trees, rooting will be done in "matchTips" in ".concon"
#    testTree <- lapply(trees, function(x) addbp(x))
  testTree <- trees
    for(j in 1:length(setBPs)){
      setBP <- setBPs[j]
      # get .concon result of specified node or all nodes of this tree,
      # as nodes are used in "refinfo" and for test trees they are just follows the nodes in refinfo
      con_results[[j]] <- lapply(tiplist2, function(x) .concon(x, testTree = testTree, setBP = setBP, reftips = reftips, from = from))
    }


  # for the processes following .concon.root2:
  for(j in 1:length(setBPs)){ # for lists of different bp results

    con_result <- con_results[[j]]
    cr <- lapply(con_result, function(xx) .conRes(xx))

    d <- as.numeric(names(cr)) # node ID
    dd <- as.numeric(unlist(cr)) # number of the same trees
    per <- dd/ddres
    data[[j]] <- data.frame(node = d, n.concor = dd, deno = ddres, percent = per)

    if(getTreeNames){
      # to get identity of trees in concordances for each node
      conTrees <- vector("list", length(con_result))
      names(conTrees) <- names(con_result)
      for(r in 1:length(conTrees)){
        t <- con_result[[r]]

        temp <- lapply(t, function(tt) tt[1,2])
        temp <- unlist(temp)

        rname <- names(t)
        w <- which(temp)
        if(length(w) > 0)
          conTrees[[r]] <- rname[w]
      }
      conTreesl[[j]] <- conTrees
    }
  } # for setBPs

  names(data) <- setBPs

  if(unresolve.ref){
    concor.temp <- list(data = data, tree = ref)
    class(concor.temp) <- "concor"
    unref <- createUnref(concor = concor.temp, bp = unbp, percent = percent)
  }

  if(getTreeNames && unresolve.ref)
    rdata <- list(refTree = ref, data = data, conTrees = conTreesl, unresolve.ref = unref)
  if(getTreeNames && !unresolve.ref)
    rdata <- list(refTree = ref, data = data, conTrees = conTreesl)
  if(!getTreeNames && unresolve.ref)
    rdata <- list(refTree = ref, data = data, unresolve.ref = unref)
  if(!getTreeNames && !unresolve.ref)
    rdata <- list(refTree = ref, data = data)

  class(rdata) <- "concor"
  return(rdata)

} # end of function



