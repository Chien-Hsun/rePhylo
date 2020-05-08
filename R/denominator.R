# identify denominators
# 2020.05.05 modified
# for correct structure for the nodes
# and moved as a separate file

###########################################################################
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param x x
#' @param ref ref
#' @param refinfo refinfo
#' @param split a character indicates the symbol used to separate
#'     species names and gene id for gene family trees. 
#'     Defaults to \code{NULL}.
#' @param c.out logical. Whether the closest outgroup species
#'     are required to exist. Defaults to \code{FALSE}. 
#' @param sub_coverage a numeric to specify the species coverage 
#'     of the clades of each node. Used as percent when <= 1, 
#'     and as count when > 1. Defaults to \code{NULL}.
#' @param max.tip an integer specifying the max number of tips
#'     for the coverage filtering. Applied only when sub_coverage 
#'     is set as percent (<= 1). Used as a cutoff for coverage as 
#'     N tips as for deep nodes having large clades.
#'     Defaults to \code{NULL}.
#' @param bp a numeric specify the minimum bp value for the node.
#'     Defaults to \code{NULL}.  
#' @export
#' @importFrom ape extract.clade
#' @importFrom phangorn Descendants
#' @seealso \code{\link{concor.node}} 
#'     \code{\link{detect_closest_out_for_GD}}
#'     \code{\link{filter_GD_by_coverage}}

# x is each single gene tree
# because ".denominator uses ref tree to do "Descendant" and "extract.clade"
# so it is no problem for trees may not be rooted

denominator <- function(x, ref, refinfo, split = NULL, 
                         c.out = FALSE, sub_coverage = NULL, 
                         max.tip = NULL, bp = NULL){
  
  allnodes <- refinfo$tdata[ ,1]
  testtips <- x$tip.label
  if(!is.null(split)){ # for gene family trees
    testtips <- lapply(testtips, function(xx) 
      unlist(strsplit(xx, split = split, fixed = TRUE))[1])
    testtips <- unlist(testtips)
  }
  
  rres <- lapply(allnodes, function(xx) 
    .deno(xx = xx, ref = ref, 
         testtree = x, testtips = testtips, 
         c.out = c.out, split = split,
         sub_coverage = sub_coverage, 
         max.tip = max.tip, bp = bp))
  res <- unlist(rres)
  names(res) <- allnodes

  return(res)
} # end of .denominator



###########################
# .deno
# main internal^2 function
#' @title internal function
#' @description Internal function of \code{rePhylo}
#' @param xx one of allnodes
#' @param ref ref
#' @param testtree testtree
#' @param testtips testtips
#' @param split a character indicates the symbol used to separate
#'     species names and gene id for gene family trees. 
#'     Defaults to \code{NULL}.
#' @param c.out logical. Whether the closest outgroup species
#'     are required to exist. Defaults to \code{FALSE}. 
#' @param sub_coverage a numeric to specify the species coverage 
#'     of the clades of each node. Used as percent when <= 1, 
#'     and as count when > 1. Defaults to \code{NULL}.
#' @param max.tip an integer specifying the max number of tips
#'     for the coverage filtering. Applied only when sub_coverage 
#'     is set as percent (<= 1). Used as a cutoff for coverage as 
#'     N tips as for deep nodes having large clades.
#'     Defaults to \code{NULL}.
#' @param bp a numeric specify the minimum bp value for the node.
#'     Defaults to \code{NULL}.  
#' @export
#' @importFrom ape extract.clade
#' @importFrom phangorn Descendants
#' @seealso \code{\link{concor.node}} 
#'     \code{\link{detect_closest_out_for_GD}}
#'     \code{\link{filter_GD_by_coverage}}
.deno <- function(xx, ref, testtree, testtips, c.out, split,
                 sub_coverage, max.tip, bp){
  nod <- xx # node of ref tip
  rlen <- length(ref$tip.label)
  
  # tip number of this current node
  # when sub_coverage is tip number (> 1)
  # correct sub_coverage if it is > n.tips 
  tmp <- ape::extract.clade(ref, node = nod)
  tmp <- tmp$tip.label
  if(!is.null(sub_coverage)){
    if(length(tmp) < sub_coverage){ 
      # it is TRUE only when sub_coverage > 1
      sub_coverage <- length(tmp)
    }
  }
  # when sub_coverage is percent (<= 1)
  # correct sub_coverage when sub_coverage > max.tip
  if(sub_coverage <= 1){
    currentN <- length(tmp)
    test <- currentN * sub_coverage
    if(!is.null(max.tip) && test > max.tip){
      sub_coverage <- max.tip
    } 
  }
  
  # identify AB tips
  des <- phangorn::Descendants(ref, node = nod, 
                               type = "children")
  ABtips <- vector("list", length(des))
  for(i in 1:length(des)){
    dd <- des[i]
    if(dd <= rlen){ # if tips
      rr <- ref$tip.label[dd]
    } # if tips
    if(dd > rlen){ # if internal nodes
      sref <- ape::extract.clade(ref, node = dd)
      rr <- sref$tip.label
    }
    ABtips[[i]] <- rr
  }
  Atips <- ABtips[[1]]
  Btips <- ABtips[[2]]
  
  # to find the subclades contains only tips of target nodes
  # with no outgroup species
  # ref tips
  if(!is.null(split)){
    tmp2 <- unlist(lapply(tmp, function(zz) 
      unlist(strsplit(zz, split = split, fixed = TRUE))[1]))
  }
  xtips <- testtips
  commontips <- which(xtips %in% tmp2)
  commontips <- testtree$tip.label[commontips]
  # now commontips is the target to see grouping in testtree
  
  if(length(commontips) == 0 | length(commontips) == 1){
    
    print(paste(xx, "FALSE", sep=" "))
    return(FALSE)

  }
  
  # check which of the mrca + all des nodes meet the requirement
  mrca <- ape::getMRCA(phy = testtree, tip = commontips)
  alldes <- phangorn::Descendants(testtree, mrca, type = "all")
  # for internal nodes but not tips
  alldes <- alldes[alldes > length(testtree$tip.label)] 
  alldes <- c(mrca, alldes) 
  # includ mrca for the analysis
  length(alldes)
  
  rm <- c()
  for(i in 1:length(alldes)){
    node <- alldes[i]
    
    split_clade <- ape::extract.clade(testtree, node)
    split_tips <- split_clade$tip.label; length(split_tips)
    split_match <- match(split_tips, commontips)
    na <- is.na(split_match)
    
    record <- FALSE
    if(all(!na)){
      # means there is no outgroup tips inside
      # then, say (A,B), see if A-B tips are in two subclade
      # that is really define the node
      des <- phangorn::Descendants(testtree, node, 
                                   type = "children")
      # alltips is to collect all tips of the entire subclade
      alltips <- vector("list",length(des))
      for(ii in 1:length(des)){
        dd <- des[ii]
        
        if(dd <= length(xtips)){ # if tips
          stips <- testtree$tip.label[dd]
          if(!is.null(split)){
            stips <- unlist(lapply(stips, function(zz) 
              unlist(strsplit(zz, split = split, fixed = TRUE))[1]))
          }
          alltips[[ii]] <- stips
          
        } # if tips
        if(dd > length(xtips)){ # if internal nodes
          sref <- ape::extract.clade(testtree, node = dd)
          stips <- sref$tip.label
          if(!is.null(split)){
            stips <- unlist(lapply(stips, function(zz) 
              unlist(strsplit(zz, split = split, fixed = TRUE))[1]))
          }
          alltips[[ii]] <- stips
        } # if internal node
      }
      
      if(any(alltips[[1]] %in% Atips) & any(alltips[[2]] %in% Btips)){
        record <- TRUE
      }
      if(any(alltips[[1]] %in% Btips) & any(alltips[[2]] %in% Atips)){
        record <- TRUE
      }
      
    } # if all grouped
    
    rm <- c(rm, record)
  } # for alldes
  subclade_id <- alldes[rm]
  length(subclade_id)
  
  if(length(subclade_id) == 0){
    print(paste(xx, "FALSE", sep=" "))
    return(FALSE)
  }
  
  # to remove nodes descendant to a valid node
  rm <- c()
  for(i in 1:length(subclade_id)){
    sn <- subclade_id[i]
    ad <- phangorn::Descendants(testtree, node = sn, 
                                type = "all")
    w <- which(subclade_id %in% ad)
    rm <- c(rm, w)
  }
  rm <- unique(rm)
  subclade_id <- subclade_id[-rm]
  length(subclade_id)
  
  if(length(subclade_id) == 0){
    print(paste(xx, "FALSE", sep=" "))
    return(FALSE)
  }
  
  # discard the subclade_id that is root
  root <- phangorn::getRoot(testtree)
  w <- subclade_id == root
  if(any(w)){
    subclade_id <- subclade_id[!w]
  }
  
  if(length(subclade_id) == 0){
    print(paste(xx, "FALSE", sep=" "))
    return(FALSE)
  }
  
  
  # to see if the subclade_node have BP > requested
  if(!is.null(bp)){
    # prepare bp table
    allnodes <- testtree$edge[testtree$edge[ ,2] > length(testtree$tip.label), 2]
    allnodes <- as.numeric(allnodes)
    allbp <- NULL
    allbp <- testtree$node.label
    if(is.null(allbp)){
      warning("bp is requested but there is no bp value in tree.")
    }
    allbp <- as.numeric(allbp[-1])
    bptab <- data.frame(nodeID = allnodes, 
                        bp = allbp)
    
    # filter subclade_node by bp
    sfun <- function(zz, bptab){
      nbp <- bptab[bptab[,1] == zz,2]
      if(is.na(nbp)){
        return(FALSE)
      } else if(nbp >= bp){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    bptest <- lapply(subclade_id, function(zz) 
      sfun(zz, bptab))
    bptest <- unlist(bptest)
    subclade_id <- subclade_id[bptest]
  }
  length(subclade_id)
  
  if(length(subclade_id) == 0){
    print(paste(xx, "FALSE", sep=" "))
    return(FALSE)
  }
  
  # for each subclade, do judgement by filter condition
  allrecord <- c()
  for(nn in 1:length(subclade_id)){
    snode <- subclade_id[nn]
    stips <- ape::extract.clade(testtree, node = snode)
    stips <- stips$tip.label
    record <- 0
    
    # filters with coverage 
    if(!is.null(sub_coverage)){
      
      if(sub_coverage <= 1){
        # valid when > % sp in each clade
        pt <- length(stips)
        pa <- length(tmp)
        per <- pt / pa
        if(per >= sub_coverage)
          record <- record +1
        
      } else if(sub_coverage > 1){
        # valid when > N. sp in each clade
        if(length(stips) >= sub_coverage)
          record <- record +1
      }
    } else {
      # did not test sub_coverage
      # i.e. each subclade has 1 tip
      record <- record +1
    }
    
    # filters with the closest out-species
    # (+1 when valid)
    # this current version ignores the relationship
    # but only ask the existence of the c.out tip(s)
    if(c.out){ 
      # identify out tips (node) in ref tree
      anc <- phangorn::Ancestors(ref, node = nod, 
                                 type = "parent");anc
      sis <- phangorn::Descendants(ref, node = anc, 
                                   type = "children")
      sis <- sis[-which(sis == nod)]
      
      # identify out tips in test tree
      sanc <- phangorn::Ancestors(testtree, node = snode, 
                                  type = "parent")
      ssis <- phangorn::Descendants(testtree, node = sanc, 
                                    type = "children");ssis
      ssis <- ssis[-which(ssis == snode)]
      if(ssis <= length(xtips)){ # if tips
        outtips <- testtree$tip.label[ssis]
      }
      if(ssis > length(xtips)){ # if internal nodes
        outtips <- ape::extract.clade(testtree, node = ssis)
        outtips <- outtips$tip.label
      }
      
      if(sis <= rlen){ # if tips
        rr <- ref$tip.label[sis]
        if(any(outtips == rr))
          record <- record +1
      }
      if(sis > rlen){ # if internal nodes
        sref <- ape::extract.clade(ref, node = sis)
        stips <- sref$tip.label
        if(any(outtips %in% stips))
          record <- record +1
      }
    } else {
      record <- record +1
    } 
    # for c.out
    
    # to summarize
    if(record == 2){
      # valid for sub_coverage +1
      # valid for c.out +1
      allrecord <- c(allrecord, "TRUE")
    } else {
      allrecord <- c(allrecord, "FALSE")
    }
  } # end of for each subclades
  
  allrecord <- as.logical(allrecord)
  if(any(allrecord)){
    res <- TRUE
  } else {
    res <- FALSE
  }
  print(paste(xx, res, sep=" "))
  return(res)
}

