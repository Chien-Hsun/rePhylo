###################################################
## filter_GD_by_coverage
###################################################
#' @title filters the GD counts by the species coverage of
#'     each of the duplicated subclades
#'
#' @description \code{filter_GD_by_coverage} is to identify for each gd whether the two 
#'     duplicated subclades contains species number more than a given level.
#'
#' @param phyto the phyto id for the analysis.
#' @param node a numeric indicating the target node id as in \code{ape}.
#'     An alternative arg of \code{phyto}.
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of gene trees to be examined. Can be unrooted.
#' @param gdtable a \code{data.frame} of the gd table file provided by Tree2GD.
#' @param sptree the species tree used to do tree reconciliation. 
#' @param sub_coverage a numeric specify the species coverage level. 
#'     If <= 1, treated as percentage. If > 1, treated as number of species.
#' @param max.tip an integer specifying the max number of tips
#'     for the coverage filtering. Applied only when sub_coverage 
#'     is set as percent (<= 1). Used as a cutoff for coverage as 
#'     N tips as for deep nodes having large clades.
#' @param split a character. The symbol used to separate species name and the 
#'     sequence number in the gene family trees.
#' @param tree_id_tab a table for the tree id during Tree2GD and the name of the 
#'     corresponding tree file. The first column is the id, and the second column
#'     is the tree names. 
#' @param up_one_node logical. Whether to choose tips of one node deeper as the 
#'     closest outgroup, as a second step of investigation. Defaults to \code{FALSE}.
#' @param plot whether to plot subclades. Defaults to \code{FALSE}.
#' @param pdfwid the width of pdf file for subtrees. Defaults to \code{10}.
#' @param pdflen the length of pdf file for subtrees. Defaults to \code{30}.
#' @param mc.cores the number of cores used for mclapply. Defaults to \code{4}.
#' @export
#' @importFrom parallel mclapply
#' @details returns a pdf plotting the gd clades, with gd pairs 
#'     in red and other species in the same clade in sptree in green.
#'     A file with prefix "Include_or_not_" indicates for each tree 
#'     and gd the requirement is fulfilled or not. A file with prefix
#'     "Trees_have_basal_tips_" shows the page number the tree fulfilled
#'     the requirement in the pdf file. 


# only for lapply
filter_GD_by_coverage <- function(
  phyto, node, trees, gdtable, sptree, split, 
  tree_id_tab, sub_coverage = NULL, max.tip = NULL,
  up_one_node = FALSE, plot = FALSE,
  pdfwid = 10, pdflen = 30, mc.cores = 2){
  
  # check if tree names matched
  t1 <- names(trees)
  t2 <- tree_id_tab[ ,2]
  m <- match(t1, t2)
  if(all(is.na(m))){
    stop("Names in tree list are different from the tree_id_tab.")
  }
  
  if(missing(phyto) && missing(node)){
    stop("Missing both phyto ids and node ids.\n")
  }
  
  if(hasArg(node)){
    phyto <- NULL
    len <- length(node)
  }
  
  if(hasArg(phyto)){
    node <- NULL
    len <- length(phyto)
  }
  
  sum <- parallel::mclapply(c(1:len), function(z) {
    .filter_coverage(phyto = phyto[z],
                     node = node[z],
                     trees = trees, 
                     gdtable = gdtable,
                     sptree = sptree,
                     sub_coverage = sub_coverage, 
                     max.tip = max.tip,
                     split = split, 
                     tree_id_tab = tree_id_tab,
                     up_one_node = up_one_node,
                     plot = plot,
                     pdfwid = pdfwid,
                     pdflen = pdflen)
  }, mc.cores = mc.cores)

  # combine all res into one table
  res_tab <- c()
  for(i in 1:length(sum)){
    rr <- sum[[i]]
    res_tab <- rbind(res_tab, rr)
  }
  row.names(res_tab) <- NULL
  
  # write result
  write.table(res_tab, "summary_sub_coverage.txt",
              quote = FALSE, row.names = FALSE,
              col.names = TRUE, sep = "\t")
  
  return(res_tab)
}



# main function but internal for lapply
#' @title Internal. Filters the GD counts by the species coverage of
#'     each of the duplicated subclades
#'
#' @description Main internal function for \code{filter_GD_by_coverage} 
#'     to identify for each gd whether the two 
#'     duplicated subclades contains species number more than a given level.
#'
#' @param phyto the phyto id for the analysis.
#' @param node a numeric indicating the target node id as in \code{ape}.
#'     An alternative arg of \code{phyto}.
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of gene trees to be examined.
#' @param gdtable a \code{data.frame} of the gd table file provided by Tree2GD.
#' @param sptree the species tree used to do tree reconciliation.
#' @param sub_coverage a numeric specify the species coverage level. 
#'     If <= 1, treated as percentage. If > 1, treated as number of species.
#' @param max.tip an integer specifying the max number of tips
#'     for the coverage filtering. Applied only when sub_coverage 
#'     is set as percent (<= 1). Used as a cutoff for coverage as 
#'     N tips as for deep nodes having large clades.
#' @param split a character. The symbol used to separate species name and the 
#'     sequence number in the gene family trees.
#' @param tree_id_tab a table for the tree id during Tree2GD and the name of the 
#'     corresponding tree file. The first column is the id, and the second column
#'     is the tree names. 
#' @param up_one_node logical. Whether to choose tips of one node deeper as the 
#'     closest outgroup, as a second step of investigation. Default \code{FALSE}.
#' @param plot whether to plot subclades. Defaults to \code{FALSE}.
#' @param pdfwid the width of pdf file for subtrees. Default \code{10}.
#' @param pdflen the length of pdf file for subtrees. Default \code{30}.
#' @export
#' @importFrom ape plot.phylo
#' @importFrom ape nodelabels
#' @importFrom ape extract.clade
#' @importFrom ape ladderize
#' @importFrom ape getMRCA
#' @importFrom phangorn Ancestors
#' @importFrom phangorn Descendants
#' @importFrom phangorn getRoot
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom graphics title
#' @keywords internal function
#' @seealso \code{\link{filter_GD_by_coverage}}

.filter_coverage<-function(
  phyto, node, trees, gdtable, sptree, sub_coverage,
  max.tip, 
  split, tree_id_tab, up_one_node = up_one_node,
  plot = plot,
  pdfwid = pdfwid, pdflen = pdflen){
  
  # prepare node id table
  phyto_tab <- make.idTable.phyto(sptree)
  
  if(!is.null(phyto) && !is.null(node)){
    cat("Both \"phyto\" and \"node\" are given. Use \"phyto\" for the analysis.\n")
  }
  
  if(!is.null(node)){
    wgdt <- phyto_tab[which(as.numeric(phyto_tab[ ,2]) == node), 1]
    wnode <- node
  }
  
  if(!is.null(phyto)){
    wgdt <- phyto
    wnode <- phyto_tab[which(phyto_tab[,1] == phyto), 2]
    wnode <- as.numeric(wnode)
  }
  
  cat(wgdt, wnode,"\n")
  
  # for one phyto_xxx
  w <- which(gdtable[ ,3] == wgdt)
  
  # if no any gd pairs
  if(length(w) == 0){
    sum<-c(wgdt,wnode,numeric(3))
    names(sum)<-c("phyto","node","meet_coverage","total","percent")
    cat(sum,sep=", ")
    cat("\n")
    
    return(sum)
  }
  
  p109 <- gdtable[w, ]
  # grep trees
  ttrees <- p109[ ,1]
  ttrees <- unique(ttrees)
  target.order <- ttrees
  m <- match(target.order, tree_id_tab[ ,1])
  target <- tree_id_tab[m,2]
  
  treenames <- names(trees)
  m <- match(target, treenames)
  ntreenames <- treenames[m]
  target.trees <- lapply(m, function(xx) 
    return(trees[[xx]]))
  names(target.trees) <- ntreenames
  
  
  # list of gd pairs
  # should find mrca of gd pairs but not all pairs in the tree
  pp <- p109[ ,2] # gd
  tt <- p109[ ,1] # tree
  dup <- duplicated(pp)
  pp <- pp[!dup]
  tt <- tt[!dup]
  plist<-vector("list",length(pp)) # gd
  names(plist)<-pp
  for(i in 1:length(pp)){
    w<-which(p109[,2] == pp[i])
    t<-lapply(w,function(xx) c(p109[xx,5],p109[xx,6]))
    t<-unlist(t)
    plist[[i]]<-t
  }
  
  # get target clade tips
  st<-ape::extract.clade(sptree,node = wnode)
  btips<-st$tip.label
  if(!is.null(sub_coverage)){
    if(length(btips) < sub_coverage){ 
      # it is TRUE only when sub_coverage > 1
      sub_coverage <- length(btips)
    }
    # when sub_coverage is percent (<= 1)
    # correct sub_coverage when sub_coverage > max.tip
    if(sub_coverage <= 1){
      currentN <- length(btips)
      test <- currentN * sub_coverage
      if(!is.null(max.tip) && test > max.tip){
        sub_coverage <- max.tip
      } 
    }
  }

  cat("All tips of the gd clade:\n")
  cat(btips,sep=", ")
  cat("\n")
  
  # to find gd pair in the target trees
  #  ntreenames
  #  target.order
  ress<-c()
  ntrees<-list()
  
  grDevices::pdf(paste("color_tip_trees_",wgdt,".pdf",sep=""),pdfwid,pdflen)
  
  for(i in 1:length(target.order)){
    options(warn = 1)
    t<-target.order[i]
    tr<-NULL
    tr<-target.trees[[i]]
    if(is.null(tr)){
      
      warning("A tree in gdtable is not present in trees object: ",t,"\n")
      
    } else {
      
      # required gd
      w<-which(tt == t)
      pairs<-lapply(w,function(xx)
        return(plist[[xx]]))
      names(pairs)<-names(plist)[w]
      
      for(y in 1:length(pairs)){
        pair<-unique(pairs[[y]])
        pn<-names(pairs)[y]
        
        # root by the most distant node
        #      options(warn=-1)
        tr <- root_dist(tr,tar = pair)
        tr <- rewrite.tree(ape::ladderize(tr,right=FALSE))
        
        options(warn = 1)
        mrca<-ape::getMRCA(tr,tip = pair);mrca
#        if(mrca != phangorn::getRoot(tr)){
        # for filter_GD_by_coverage,
        # mrca == root is fine.
          des<-phangorn::Descendants(tr,node = mrca,
                                     type = "children");des
          # filter for the two des 
          # == the two duplicated subclades
          tmp <- 0
          for(ndes in 1:length(des)){
            sis <- des[ndes]
            
            if(sis <= length(tr$tip.label)){
              # tips
              stips<-tr$tip.label[sis]
            } else { # nodes
              str<-ape::extract.clade(tr,node = sis)
              stips<-str$tip.label
            }
            
            # see the coverage of subclade sp
            if(!is.null(split)){
              stips<-lapply(stips,function(xx)
                unlist(strsplit(xx,split = split,fixed = TRUE))[1])
              stips<-unlist(stips)
            }
            check<-stips %in% btips
            if(sub_coverage <= 1){ # percent
              
              blen <- length(btips)
              len <- length(which(check))
              len <- len / blen
              if(len >= sub_coverage){
                tmp <- tmp + 1
              }
              
            } else if(sub_coverage > 1){ # count
              
              len <- length(which(check))
              if(len >= sub_coverage){ # number of covered tips
                tmp <- tmp + 1
              }
            }
          } # for the 2 des
          
          # summarize result for the two nodes
          # tmp == 2 means both subclade meet the requirement
          if(tmp == 2){
            res<-c(ntreenames[i], pn, "TRUE")
          } else {
            res<-c(ntreenames[i], pn, "FALSE")
          }
          ress<-rbind(ress,res)
          
          # plot trees
          # in "detect_out" extracts subtree include up-two-nodes
          # here can just extract subtree
          atree<-ape::extract.clade(phy = tr,node = mrca)
          # prepare tip colors
          tipcol<-rep("black",length(atree$tip.label))
          tips<-atree$tip.label
          if(!is.null(split)){
            tips<-lapply(tips,function(xx)
              unlist(strsplit(xx,split = split,fixed = TRUE))[1])
            tips<-unlist(tips)
          }
          # color basal tips in green
          w<-which(tips %in% btips);w
          if(length(w) > 0){tipcol[w]<-"green"}
          # color pairs in red
          w<-which(atree$tip.label %in% pair);w
          if(length(w) > 0){tipcol[w]<-"red"}
          # plot tree
          ape::plot.phylo(atree,tip.color = tipcol,
                          use.edge.length = FALSE,
                          node.depth = 2)
          ape::nodelabels(text = atree$node.label,
                          bg = "white")
          graphics::title(paste("tree = ", ntreenames[i],
                                ", gd = ", pn, sep=""))
          atree$tip.col <- tipcol
          alist <- list(atree)
          names(alist) <- paste(ntreenames[i], ", gd = ", 
                              pn, sep="")
          ntrees <- append(ntrees,alist)
#        } else {
#          war <- paste(ntreenames[i], ", gd ", 
#                       pn, " has mrca == root.", sep="")
#          warning(war)
#          res <- c(war, " ", " ")
#          ress <- rbind(ress, res)
#          trr <- list(tr)
#          names(trr) <- paste(ntreenames[i], ", gd = ", 
#                              pn, sep="")
#          ntrees <- append(ntrees, trr)
#        }
      } # for gd
      
    } # if tree is not NULL

  } # for trees
  grDevices::dev.off()
  #
  saveRDS(ntrees,paste("ntrees_",wgdt,".rds",sep=""))
  write.table(ress,paste("Meet_coverage_or_not_",wgdt,".txt",sep=""),
              quote = FALSE,col.names = TRUE,row.names = FALSE)
  # ress contains, for each gd, TRUE or FALSE to have the outgroup
  # the first column is the tree id, the second column is TREU/FALSE
  # thus the first column may be duplicated
  
  tar<-ress[,3]
  w<-which(tar == "TRUE")
  lw<-length(w) # length of having closest outgroups
  lt<-length(tar) # length of all gd
  sum<-c(wgdt,wnode,lw,lt,round(lw/lt,3))
  names(sum)<-c("phyto","node","meet_coverage","total","percent")
  cat(sum,sep=", ")
  cat("\n")
  r<-ress[w,1]
  
  # write a table
  page<-w
  trn<-r
  dd<-data.frame(page=page,tree=trn)
  write.table(dd,paste("Trees_meet_coverage_",wgdt,".txt",sep=""),
              sep="\t",quote = FALSE,col.names = TRUE,
              row.names = FALSE)
  
  return(sum)
}
