##################################################################################
# internal functions
##################################################################################
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param tree tree
#' @param subtips subtips
#' @export
#' @noRd
#' @keywords internal
.compare.edge.len <- function(tree, subtips){
  bls <- numeric(length(subtips))
  for(j in 1:length(subtips)){
    stt <- subtips[j]
    tt <- tree$tip.label
    ww <- which(tt == stt)
    bls[j] <- tree$edge.length[tree$edge[ ,2] == ww]
  }
  return(bls)
}


##################################################################################
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param zz zz
#' @param n n
#' @param all.dedup all.dedup
#' @export
#' @noRd
#' @keywords internal
.sum.num <- function(zz, n, all.dedup){
  cons <- zz
  w <- which(all.dedup[ ,n] == zz)
  len <- length(w)
  return(len)
}


##################################################################################
# check all tips of all trees, with ref$tip.label or tips in taxa table
# if there are tips lacking in all test trees comparing to those in ref/taxa table
# or vise versa
# then will remove these tips (which will cause trouble in the main codes)
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param trees trees
#' @param ref ref
#' @param taxa taxa
#' @export
#' @noRd
#' @importFrom ape drop.tip
#' @keywords internal
.match.tiplabels <- function(trees, ref = NULL, taxa = NULL){

  allt <- lapply(trees, function(x) return(x$tip.label))
  allt <- unlist(allt)
  allt <- as.character(allt[!duplicated(allt)])

  if(!is.null(ref)){
    # for ref, all tips in ref should be present in allt (all tips from test trees)
    # and tips in all test trees should also present in ref (not having extra tips)
    rtips <- ref$tip.label
    if(!all(allt %in% rtips)){
      warning("There are extra tips in test trees not present in ref tree. Will remove these tips in all test trees.")
      m <- match(allt, rtips)
      lacktips <- allt[is.na(m)]
      tname<-names(trees)
      trees <- lapply(trees, function(x) ape::drop.tip(phy = x, tip = lacktips))
      names(trees)<-tname
    }
    if(!all(rtips %in% allt)){
      warning("There are extra tips in ref tree not present in any of the test trees. Extra tips will be removed in ref tree.")
      m <- match(rtips, allt)
      lacktips <- rtips[is.na(m)]
      ref <- ape::drop.tip(phy = ref, tip = lacktips)
    }
    res <- list(trees = trees, ref = ref)
    return(res)
  }

#  if(is.null(ref) && !is.null(taxa)){
  # maybe better not to put priority for ref over taxa table in this function
  # because it will be used in more than one function
  if(!is.null(taxa)){
    taxat <- as.character(taxa[ ,1])

    if(!all(taxat %in% allt)){
      # there might be tips in trees not included in taxa (constraints), it is normal
      # but if there are tips in taxa not in allt (all tips of all trees), this is a problem
      warning("There is tip(s) in taxa table not shown in all the test trees. Will remove lacking tips in taxa table.")
      m <- match(taxat, allt)
      taxa <- taxa[!is.na(m), ]
    }
    res <- list(trees = trees, taxa = taxa)
    return(res)
  }
}




#########################################################################
# internal function
# finds mrca node of the specified constraint for a test tree
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param constr constr
#' @param tree tree
#' @param taxa taxa
#' @param level level
#' @param bp a logical indicating whether to retain
#'     the original bp values in the input tree.
#'     Currently desinged for backboneBP.
#' @export
#' @noRd
#' @importFrom ape getMRCA
#' @keywords internal
#' @seealso \code{\link{backboneBP}}

.getmnode <- function(constr, tree, taxa, level, bp = FALSE){

  ww <- which(taxa[ ,level] == constr);length(ww)
  wtips <- taxa[ww,1];length(wtips)
  wtips <- as.character(wtips)
  yy <- which(wtips %in% tree$tip.label)

  if(length(yy) == 0){

    res <- NULL

  } else {

    wtips <- wtips[yy]; length(wtips) # remove those not in the trees (after drop.tip)
    # removing tips is mainly required by backboneBP ??

    tree <- root.dist(z = tree, tar = wtips, bp = bp)
    # because wtips should be characters, so can use directly
    mrca <- ape::getMRCA(tree, tip = wtips); mrca
    res <- list(tree = tree, mrca = mrca)

  }

  return(res)
}




#########################################################################
# internal function
# finds mrca node of the specified constraint for a test tree
#' @title internal function
#' @description Identify backbone nodes based on \code{taxa} and
#'     make other bp values as NA. Internal function of \code{rePhylo}.
#' @param tre a tree
#' @param mrcalist mrcalist
#' @export
#' @importFrom phangorn Descendants
#' @keywords internal
#' @seealso \code{\link{backboneBP}}
# prepare function to make the bp of node not backbone n as NA
# input x as testSet , the number of list in length of trees
# the idea is to root the tree that have previously edited bp
# by the root based on mrcalist
# and do cycle for the same object to remove
# all the bp that are not backbone nodes
select.backbone <- function(tre, mrcalist){

  if(!is.null(mrcalist)){

    mrca <- mrcalist$mrca
    tr <- mrcalist$tree # from .getmnode
#    rr <- phangorn::getRoot(tr)
#    dd <- phangorn::Descendants(tr, node = rr, type = "children")
#    dd <- which(dd <= length(tr$tip.label))
#    dd <- tr$tip.label[dd]
#    tr <- ape::root.phylo(phy = tre, outgroup = dd,
#                          resolve.root = TRUE, edgelabel = TRUE)
#    tr <- rewrite.tree(tr)

    if(!is.null(mrca)){
      des <- phangorn::Descendants(tr, mrca, type = "all");des
      ee <- tr$edge[tr$edge[, 2] > length(tr$tip.label), 2]
      wi <- which(ee %in% des); length(wi)
      # because root has a node.label but absent in edge matrix
      # so need to +1
      wi <- wi+1
      nl <- tr$node.label
      nl[wi] <- NA
      # the below two lines is to replace "addbp"
      # that make node labels of "Root" or "" (not numbers) as NA
      nl <- as.numeric(nl)
      nl <- as.character(nl)
      tr$node.label <- nl

    } # end of if mrca != NULL
    # else {
    # if test.mrca == NULL , indicating there might be missing tips in trees
    # and then there is only one tip of the constraint in the tree ??
    # in this case there is no node value to be set as NA , so just do nothing
    #     }
  } else {
    tr<-tre
    cat("There is NULL for mrca")
  }

  return(tr)
} # end of select.backbone






# internal functions for "concor.node"
# new version: judge based on dres (denominotor pools)

#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param z z
#' @param node node
#' @importFrom ape extract.clade
#' @export
#' @noRd
#' @seealso \code{\link{concor.node}}

# using "testTree" as z's list, so z is each tree
# because ".getInfo" uses ref tree to do "extract.clade"
# so it is no problem
.getInfo <- function(z, node){ # to get node ids, bps, tips for each node

  allnodes <- z$edge[z$edge[ ,2] > length(z$tip.label),2]

  if(!is.null(node))
    nodes <- node
  else
    nodes <- allnodes

  # 1st list: table with nodeID+bp values
  mat <- matrix(nrow = length(allnodes), ncol = 2)
  colnames(mat) <- c("nodeID", "bp")
  mat[ ,1] <- as.numeric(allnodes)
  bp <- z$node.label
  if(is.null(bp)){
    mat[ ,2] <- NA
  } else {
    bp <- bp[-1]
    mat[ ,2] <- as.numeric(bp)
  }

  mat <- mat[mat[ ,1] %in% nodes, ]

  if(inherits(mat, "numeric"))
    mat <- t(as.matrix(mat))
  colnames(mat) <- c("nodeID", "bp")

  # 2nd list: all tips for each node
  alist <- as.list(nodes)
  names(alist) <- nodes
  # to get tips for all nodes
  for(q in 1:length(nodes)){
    no <- nodes[q]

    stree <- ape::extract.clade(z, node = no)
    stips <- stree$tip.label

    alist[[q]] <- stips
  }
  tlist <- list(tdata = mat, stips = alist)

  return(tlist)
}




###########################################################################
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param x x
#' @param ref ref
#' @param refinfo refinfo
#' @export
#' @noRd
#' @importFrom ape extract.clade
#' @importFrom phangorn Descendants
#' @seealso \code{\link{concor.node}}

# x is each single gene tree
# because ".denominator uses ref tree to do "Descendant" and "extract.clade"
# so it is no problem
.denominator <- function(x, ref, refinfo){

  allnodes <- refinfo$tdata[ ,1]
  testtips <- x$tip.label
  rlen <- length(ref$tip.label)

  res <- as.list(allnodes)
  names(res) <- allnodes

  deno <- function(xx, ref, testtips){
    nod <- xx
    des <- phangorn::Descendants(ref, node = nod, type = "children");des
    record <- 0
    for(d in 1:length(des)){
      dd <- des[d]
      if(dd <= rlen){ # if tips
        rr <- ref$tip.label[dd]
        if(any(testtips == rr))
          record <- c(record, 1)
      }
      if(dd > rlen){ # if internal nodes
        sref <- ape::extract.clade(ref, node = dd)
        stips <- sref$tip.label
        if(any(stips %in% testtips))
          record <- c(record, 1)
      }
    }
    sum <- sum(record)
    if(sum >= 2)
      res <- TRUE
    else res <- FALSE

    return(res)
  }

  rres <- lapply(allnodes, function(xx) deno(xx, ref, testtips))
  res <- unlist(rres)
  names(res) <- allnodes
  return(res)
} # end of .denominator




# z is tiplist's N object, equals to tips of each node in original ref tree
# results are for all or specified nodes in a tree
# INPUTS: one node on ref, all testTrees

# for rooted trees (no rooting process within this package)
# concon.ori <- function(x, testTree, setBP, ref){
#  refTips <- x$tips

# input are: testTree as a list of many trees, refTips is tips of one node
#  matchTips_res <- lapply(testTree, function(xx) .matchTips(xx, refTips, setBP, ref))

#  return(matchTips_res)
#}

# for two steps rooting of concordance

# .concon should be the same for either rooted or unrooted tree
# because this is after considering denominators, and this issue occurrs due the nature of single gene trees
# that some tips might lost in single gene trees
# so this is unrelated to the rooting process within the codes
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param x x
#' @param testTree testTree
#' @param setBP setBP
#' @param reftips reftips
#' @param from from
#' @export
#' @noRd
#' @seealso \code{\link{concor.node}}, \code{\link{backboneBP}}, \code{\link{concor.trees}}

.concon <- function(x, testTree, setBP, reftips, from){
  refTips <- x$tips
  trs <- x$trees.w.node

  if(from == "concor.trees"){
    matchTips_res <- lapply(testTree, function(xx)
      .matchTips(zz = xx, refTips = refTips, setBP = setBP, reftips = reftips))
    return(matchTips_res)
  }

  if(from == "concor.node"){
    trs <- x$trees.w.node

    num <- vector("list", length(testTree))
    names(num) <- names(testTree)
    for(i in 1:length(num)){
      w <- which(trs == i)
      if(length(w) == 0)
        w <- NULL
      num[[i]] <- list(after = w, ori = i)
    }

    trs <- as.list(trs)
    trs <- lapply(trs, function(x) return(testTree[[x]]))
    # input are: testTree as a list of many trees, refTips is tips of one node
    matchTips_res <- lapply(trs, function(xx) 
      .matchTips(zz = xx, refTips = refTips, setBP = setBP, 
                 reftips = reftips))

    ttx <- function(xx, matchTips_res){
      after <- xx$after
      if(!is.null(after))
        rr <- matchTips_res[[after]]
      else{
        record <- tn <- FALSE
        rr <- data.frame(tn, record)
      }
      return(rr)
    }
    res <- lapply(num, function(xxx) ttx(xxx, matchTips_res))
    return(res)
  }
}


# input as z: one testTree, and refTips: the tips of one node of ref tree
# INPUTS: one node on ref, one testTree
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param zz a tree (class phylo)
#' @param refTips target tips to be matched
#' @param setBP setBP
#' @param reftips all tips in ref tree
#' @export
#' @noRd
#' @importFrom ape getMRCA
#' @importFrom ape extract.clade
#' @importFrom phangorn getRoot
#' @seealso \code{\link{concor.node}}, \code{\link{backboneBP}}, \code{\link{concor.trees}}

# 2019.05.18
# has implemented "root.dist" internal function to root test tree
.matchTips <- function(zz, refTips, setBP, reftips){

  # first root the test tree based on refTips
  isroot<-is.rooted(zz)
  if(!isroot){
    zz <- root.dist(z = zz, tar = refTips)
  }

  if(is.null(zz)){

    record <- FALSE
    tn <- FALSE

  } else {

    # to remove tips (in ref tree, which contain the msot tips) not included in test trees
    alltips <- zz$tip.label
    #  reftips <- ref$tip.label
    match <- match(reftips, alltips)
    lacktips <- reftips[is.na(match)]

    which <- which(refTips %in% lacktips)
    if(length(which) > 0){
      refTips <- refTips[-which]
    }

    mn <- ape::getMRCA(phy = zz, tip = refTips)
    if(!is.null(mn)){
      stree <- ape::extract.clade(phy = zz, node = mn)
      mrcatips <- stree$tip.label

      record <- FALSE
      tn <- FALSE
      if(all(mrcatips %in% refTips) && all(refTips %in% mrcatips)){
        record <- TRUE
        tn <- mn
      }
    } else {
      # if there is only one tip in getMRCA, NULL is returned
      # NOT record because we did not count for tips (only see internal nodes)
      record <- FALSE
      tn <- FALSE
    }
  }



  if(record){ # to judge bp value

    root <- phangorn::getRoot(zz)
    if(root != tn){
      bps <- as.numeric(zz$node.label[-1])
      tm <- as.numeric(zz$edge[zz$edge[ ,2] > length(zz$tip.label),2])
      ww <- which(tm == tn)
      bp <- NULL
      if(!is.null(ww) & length(ww) != 0)
        bp <- bps[ww]
      if(!is.null(bp) && !is.na(bp)){
        # so if there is empty bp, then will record (recognize)
        # this node regardless of bp values
        bp <- as.numeric(bp)
        if(bp < setBP){record <- FALSE}
      }
    }
  }

  res <- data.frame(tn, record)
  # we have the node ID if there is matched tips, and c(0, 0) if there is no any matched clades
  # here "tn" is the node in the tree that rooted by the most distant tip
  # using "root.dist" internal function
  # so will be different in different test trees

  return(res)
}





#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param z z
#' @export
#' @noRd
#' @keywords internal
#' @seealso \code{\link{concor.node}}, \code{\link{backboneBP}}, \code{\link{concor.trees}}

.conRes <- function(z){

  if(!is.null(z)){
    rr <- 0
    for(i in 1:length(z)){
      tar <- z[[i]]
      if(is.null(tar)){rr <- 0}
      else{
        tt <- tar[1,2]
        if(tt == TRUE) # w is the ones that are concordant
          # "TRUE" become as characters when using rbind above, so == "TRUE" is required
          rr <- c(rr, 1)
      }
    }
    ww <- sum(rr)

    return(ww)

  } else {return(0)}
}





# 2019.05.18 new internal rooting function
# to root tree using the most distant tip as temporary root

#' @title internal function
#' @description Re-numbering the tips and nodes of ape tree objects.
#'     Internal function of \code{rePhylo}.
#' @param z an object of class "\code{phylo}"
#' @param tar target tips in the current analysis
#' @param bp a logical indicating whether to retain
#'     the original bp values in the input tree.
#'     Currently desinged for backboneBP.
#' @importFrom ape dist.nodes
#' @importFrom ape root.phylo
#' @importFrom phangorn Descendants
#' @importFrom phangorn getRoot
#' @export
#' @keywords internal function
#' @seealso \code{\link{cladeFilter}}, \code{\link{concor.node}}, \code{\link{backboneBP}}, \code{\link{concor.trees}}

root.dist<-function(z, tar, bp = FALSE){
  # z is the tree (one tree)
  # tar is the target in the "current" analysis (one point, one node)
  rt<-z

  if(inherits(tar,"character"))
    tartip<-which(rt$tip.label %in% tar)
  if(inherits(tar,"numeric"))
    tartip<-tar
  # now "tartip" is the number id of target tips

#  dd<-ape::cophenetic.phylo(x = rt)
  dd<-ape::dist.nodes(x = rt)
  dd2<-NULL
  if(length(tartip) == length(rt$tip.label)){
    # indicates all tips are as one group
    # don't need to do rooting
    return(rt)
  }

  if(length(tartip) == (length(rt$tip.label)-1)){
    # there is only one tip not in the group
    # so just root with this tip
    ww2<-c(1:length(rt$tip.label))[-tartip]
    dd2<-dd[ww2,]
    dd2<-data.frame(dd2)
  } else if(length(tartip) > 1) {
    tt<-rt$tip.label[tartip]
    dd2<-dd[,tartip];ncol(dd2);nrow(dd2)
    dd2<-dd2[-tartip,];ncol(dd2);nrow(dd2)

    max2<-apply(dd2,2,max);max2
    w2<-lapply(c(1:ncol(dd2)),function(x) which(dd2[,x] == max2[x]));w2
    ww2<-as.numeric(unlist(lapply(w2,names)));ww2
    ww2<-ww2[!duplicated(ww2)][1]
  } else if(length(tartip) == 1) {
    tt<-rt$tip.label[tartip]
    dd2<-dd[-tartip,]
    dd2<-dd2[,tartip]

    w2<-which(dd2 == max(dd2));w2
    #    w2<-lapply(c(1:ncol(dd2)),function(x) which(dd2[,x] == max2[x]));w2
    ww2<-as.numeric(names(w2))[1];ww2
    #    ww2<-ww2[!duplicated(ww2)]
  }

  if(is.null(dd2)){
    # in this test tree there is no tips belongs to this target constraint
    return(NULL)
  }

  if(ww2 > length(rt$tip.label)){ # is an internal node (from dist.nodes)
    rt2<-ape::extract.clade(phy = rt, node = ww2)
    root<-rt2$tip.label
  } else {
    if(length(ww2) == 1){
      root<-rt$tip.label[ww2]
    } else {
      dd3<-dd2[ww2,]
      w3<-apply(dd3,1,sum)
      ww3<-which(w3 == max(w3))
      ww3<-as.numeric(names(ww3))
      root<-rt$tip.label[ww3]
    }
  }

  rr <- phangorn::getRoot(rt)
  rr <- phangorn::Descendants(x = rt, node = rr, type = "children")
  rr <- rr[rr <= length(rt$tip.label)]
  if(length(rr) > 0){
    if(all(rt$tip.label[rr] %in% root) & bp){
      return(rt)
    }
  }
  rrt <- ape::root.phylo(phy = rt, outgroup = root,
                         resolve.root = TRUE, edgelabel = TRUE)
  options(warn=-1)
  rrt$node.label <- as.character(as.numeric(rrt$node.label,warn=FALSE))
  rrt <- rewrite.tree(rrt)
  options(warn=1)
  
  return(rrt)
}


