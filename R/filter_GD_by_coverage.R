###################################################
## filter_GD_by_coverage
###################################################
#' @title filters the GD counts by the species coverage of
#'     each of the duplicated subclades
#'
#' @description \code{filter_GD_by_coverage} is to identify for each gd whether the two 
#'     duplicated subclades contains species number more than a given level.
#'
#' @param phyto_node a table with the first column as the phyto id given by Tree2GD and 
#'     the second column as the node id from \code{ape}.
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of gene trees to be examined. Can be unrooted.
#' @param gdtable a \code{data.frame} of the gd table file provided by Tree2GD.
#' @param sptree the species tree used to do tree reconciliation. 
#' @param sub_coverage a numeric specify the species coverage level. 
#'     If <= 1, treated as percentage. If > 1, treated as number of species.
#' @param split a character. The symbol used to separate species name and the 
#'     sequence number in the gene family trees.
#' @param tree_id_tab a table for the tree id during Tree2GD and the name of the 
#'     corresponding tree file. The first column is the id, and the second column
#'     is the tree names. 
#' @param up_one_node logical. Whether to choose tips of one node deeper as the 
#'     closest outgroup, as a second step of investigation. Default \code{FALSE}.
#' @param pdfwid the width of pdf file for subtrees. Default \code{10}.
#' @param pdflen the length of pdf file for subtrees. Default \code{30}.
#' @export
#' @details returns a pdf plotting the gd clades, with gd pairs 
#'     in red and other species in the same clade in sptree in green.
#'     A file with prefix "Include_or_not_" indicates for each tree 
#'     and gd the requirement is fulfilled or not. A file with prefix
#'     "Trees_have_basal_tips_" shows the page number the tree fulfilled
#'     the requirement in the pdf file. 


# only for lapply
filter_GD_by_coverage <- function(
  phyto_node, trees, gdtable, sptree, split, 
  tree_id_tab, sub_coverage = 2,
  up_one_node = FALSE, 
  pdfwid = 10, pdflen = 30){
  
  # check if tree names matched
  t1 <- names(trees)
  t2 <- tree_id_tab[ ,2]
  m <- match(t1, t2)
  if(all(is.na(m))){
    stop("Names in tree list are different from the tree_id_tab.")
  }
  
  sum<-lapply(1:nrow(phyto_node),function(z) 
    .filter_coverage(x = z, 
                phyto_node = phyto_node,
                trees = trees, 
                gdtable = gdtable,
                sptree = sptree,
                sub_coverage = sub_coverage, 
                split = split, 
                tree_id_tab = tree_id_tab,
                up_one_node = up_one_node,
                pdfwid = pdfwid,
                pdflen = pdflen))

  return(sum)
}



# main function but internal for lapply
#' @title Internal. Filters the GD counts by the species coverage of
#'     each of the duplicated subclades
#'
#' @description Main internal function for \code{filter_GD_by_coverage} 
#'     to identify for each gd whether the two 
#'     duplicated subclades contains species number more than a given level.
#'
#' @param x the row of phyto_node to be tested
#' @param phyto_node a table with the first column as the phyto id given by Tree2GD and 
#'     the second column as the node id from \code{ape}.
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of gene trees to be examined.
#' @param gdtable a \code{data.frame} of the gd table file provided by Tree2GD.
#' @param sptree the species tree used to do tree reconciliation.
#' @param sub_coverage a numeric specify the species coverage level. 
#'     If <= 1, treated as percentage. If > 1, treated as number of species.
#' @param split a character. The symbol used to separate species name and the 
#'     sequence number in the gene family trees.
#' @param tree_id_tab a table for the tree id during Tree2GD and the name of the 
#'     corresponding tree file. The first column is the id, and the second column
#'     is the tree names. 
#' @param up_one_node logical. Whether to choose tips of one node deeper as the 
#'     closest outgroup, as a second step of investigation. Default \code{FALSE}.
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
  x, phyto_node, trees, gdtable, sptree, sub_coverage,
  split, tree_id_tab, up_one_node = up_one_node,
  pdfwid = pdfwid, pdflen = pdflen){
  
  wgdt <- phyto_node[x,1]
  # for one phyto_xxx
  w <- which(gdtable[ ,3] == wgdt)
  length(w)
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
  w<-which(phyto_node[,1] == wgdt)
  des<-as.numeric(phyto_node[w,2])
  st<-ape::extract.clade(sptree,node = des)
  btips<-st$tip.label
  if(length(btips) < sub_coverage){
    sub_coverage <- length(btips)
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
        tr<-root.dist(tr,tar = pair)
        tr<-rewrite.tree(ape::ladderize(tr,right=FALSE))
        
        options(warn = 1)
        mrca<-ape::getMRCA(tr,tip = pair);mrca
        if(mrca != phangorn::getRoot(tr)){
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
            stips<-lapply(stips,function(xx)
              unlist(strsplit(xx,split = split,fixed = TRUE))[1])
            stips<-unlist(stips)
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
          tips<-lapply(tips,function(xx)
            unlist(strsplit(xx,split = split,fixed = TRUE))[1])
          tips<-unlist(tips)
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
        } else {
          war <- paste(ntreenames[i], ", gd ", 
                       pn, " has mrca == root.", sep="")
          warning(war)
          res <- c(war, " ", " ")
          ress <- rbind(ress, res)
          trr <- list(tr)
          names(trr) <- paste(ntreenames[i], ", gd = ", 
                              pn, sep="")
          ntrees <- append(ntrees, trr)
        }
      } # for gd
      
    } # if tree is not NULL

  } # for trees
  grDevices::dev.off()
  #
  saveRDS(ntrees,paste("ntrees_",wgdt,".rds",sep=""))
  write.table(ress,paste("Include_or_not_",wgdt,".txt",sep=""),
              quote = FALSE,col.names = TRUE,row.names = FALSE)
  # ress contains, for each gd, TRUE or FALSE to have the outgroup
  # the first column is the tree id, the second column is TREU/FALSE
  # thus the first column may be duplicated
  
  tar<-ress[,3]
  w<-which(tar == "TRUE")
  lw<-length(w) # length of having closest outgroups
  lt<-length(tar) # length of all gd
  sum<-c(wgdt,lw,lt,round(lw/lt,3))
  names(sum)<-c("wgd_node","having_out","total","percent")
  cat(sum,sep=", ")
  cat("\n")
  r<-ress[w,1]
  
  # write a table
  page<-w
  trn<-r
  dd<-data.frame(page=page,tree=trn)
  write.table(dd,paste("Trees_have_basal_tips_",wgdt,".txt",sep=""),
              sep="\t",quote = FALSE,col.names = TRUE,
              row.names = FALSE)
  
  return(sum)
}
