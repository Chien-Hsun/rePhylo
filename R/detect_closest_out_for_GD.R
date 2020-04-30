###################################################
## detect_closest_out_for_GD
###################################################
#' @title Detect the existence of the closest outgroups for each gd
#'
#' @description \code{detect_closest_out_for_GD} is to find out whether each gd 
#'     has tips of its sister clade as the closest outgroup, in order to verify 
#'     the support of gd. The returning object contains a simple summary for the 
#'     number and percent of gd that has the closest outgroup. This function will 
#'     also produce pdf files plotting subclade of each gd, with tips colored, 
#'     for easier manual investigation when required.
#'
#' @param phyto_node a table with the first column as the phyto id given by Tree2GD and 
#'     the second column as the node id from \code{ape}.
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of gene trees to be examined. Can be unrooted.
#' @param gdtable a \code{data.frame} of the gd table file provided by Tree2GD.
#' @param sptree the species tree used to do tree reconciliation. 
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
#' @details For each gd, detect whether there are closest outgroups
#'    exactly next to the duplicated node, that can actually define 
#'    the gd at the indicated node. Otherwise, gd will be possible 
#'    to be at other (depper) nodes as well.
#' @examples
#' \dontrun{
#' res<-
#' detect_closest_out_of_GD(
#' phyto_node,trees = trees,gdtable = gdtable,
#' sptree = sptree,split=".",tree_id_tab = tree_id_tab,
#' up_one_node=FALSE,pdfwid = 20,pdflen = 100)
#' 
#' res
#' }

# only for lapply
detect_closest_out_for_GD <- function(
  phyto_node, trees, gdtable, sptree, split, 
  tree_id_tab, up_one_node = FALSE, 
  pdfwid = 10, pdflen = 30){
  
  # check if tree names matched
  t1 <- names(trees)
  t2 <- tree_id_tab[ ,2]
  m <- match(t1, t2)
  if(all(is.na(m))){
    stop("Names in tree list are different from the tree_id_tab.")
  }
  
  sum<-lapply(1:nrow(phyto_node),function(z) 
    .detect_out(z, phyto_node,trees, 
                gdtable,sptree,
                split,
                tree_id_tab,
                up_one_node = up_one_node,
                pdfwid = pdfwid,
                pdflen = pdflen))

  return(sum)
}



# main function but internal for lapply
#' @title Internal. Detect the existence of the closest outgroups for each gd
#'
#' @description Main internal function for \code{detect_closest_out_for_GD} 
#'     that is to find out whether each gd 
#'     has tips of its sister clade as the closest outgroup, in order to verify 
#'     the support of gd. The returning object contains a simple summary for the 
#'     number and percent of gd that has the closest outgroup. This function will 
#'     also produce pdf files plotting subclade of each gd, with tips colored, 
#'     for easier manual investigation when required.
#'
#' @param x the row of phyto_node to be tested
#' @param phyto_node a table with the first column as the phyto id given by Tree2GD and 
#'     the second column as the node id from \code{ape}.
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of gene trees to be examined.
#' @param gdtable a \code{data.frame} of the gd table file provided by Tree2GD.
#' @param sptree the species tree used to do tree reconciliation.
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
#' @details For each gd, detect whether there are closest outgroups
#'    exactly next to the duplicated node, that can actually define 
#'    the gd at the indicated node. Otherwise, gd will be possible 
#'    to be at other (depper) nodes as well.
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
#' @seealso \code{\link{detect_closest_out_for_GD}}

.detect_out<-function(
  x, phyto_node, trees, gdtable, sptree,
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
  
  
  # get basal tips
  w<-which(phyto_node[,1] == wgdt)
  bnode<-as.numeric(phyto_node[w,2])
  anc<-phangorn::Ancestors(sptree,node = bnode,type = "parent")
  des<-phangorn::Descendants(sptree,node = anc,type = "children")
  des<-des[-which(des == bnode)]
  if(des <= length(sptree$tip.label)){
    btips<-sptree$tip.label[des]
  } else {
    st<-ape::extract.clade(sptree,node = des)
    btips<-st$tip.label
  }
  if(up_one_node){
    bnode<-anc
    anc<-phangorn::Ancestors(sptree,node = bnode,type = "parent")
    des<-phangorn::Descendants(sptree,node = anc,type = "children")
    des<-des[-which(des == bnode)]
    if(des <= length(sptree$tip.label)){
      btips<-sptree$tip.label[des]
    } else {
      st<-ape::extract.clade(sptree,node = des)
      btips<-st$tip.label
    }
  }
  cat("basal_tips:\n")
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
        
        options(warn=1)
        mrca<-ape::getMRCA(tr,tip = pair);mrca
        if(mrca != phangorn::getRoot(tr)){
          anc<-phangorn::Ancestors(tr,node = mrca,type = "parent");anc
          des<-phangorn::Descendants(tr,node = anc,type = "children");des
          sis<-des[-which(des == mrca)];sis
          
          if(sis <= length(tr$tip.label)){ # tips
            stips<-tr$tip.label[sis]
          } else { # nodes
            str<-ape::extract.clade(tr,node = sis)
            stips<-str$tip.label
          }
          # see if contain basal tips
          stips<-lapply(stips,function(xx)
            unlist(strsplit(xx,split = split,fixed = TRUE))[1])
          stips<-unlist(stips)
          check<-stips %in% btips
          if(any(check)){ # having basal tips
            res<-c(ntreenames[i],"TRUE")
          } else {
            res<-c(ntreenames[i],"FALSE")
          }
          ress<-rbind(ress,res)
          
          # plot trees
          # to extract subtree include up-two-nodes
          if(anc != phangorn::getRoot(tr))
            anc<-phangorn::Ancestors(tr,node = anc,type = "parent")
          atree<-ape::extract.clade(phy = tr,node = anc)
          # color pairs in red
          tipcol<-rep("black",length(atree$tip.label))
          w<-which(atree$tip.label %in% pair);w
          if(length(w) > 0){tipcol[w]<-"red"}
          tips<-atree$tip.label
          tips<-lapply(tips,function(xx)
            unlist(strsplit(xx,split = split,fixed = TRUE))[1])
          tips<-unlist(tips)
          # color basal tips in green
          w<-which(tips %in% btips);w
          if(length(w) > 0){tipcol[w]<-"green"}
          # plot tree
          ape::plot.phylo(atree,tip.color = tipcol,use.edge.length = FALSE,
                          node.depth = 2)
          ape::nodelabels(text = atree$node.label,bg = "white")
          graphics::title(paste("tree = ",ntreenames[i],", gd = ",pn,sep=""))
          atree$tip.col<-tipcol
          alist<-list(atree)
          names(alist)<-paste(ntreenames[i],", gd = ",pn,sep="")
          ntrees<-append(ntrees,alist)
        } else {
          war<-paste(ntreenames[i],", gd ",pn," has mrca == root.",sep="")
          warning(war)
          res<-c(war," ")
          ress<-rbind(ress,res)
          trr<-list(tr)
          names(trr)<-paste(ntreenames[i],", gd = ",pn,sep="")
          ntrees<-append(ntrees,trr)
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
  
  tar<-ress[,2]
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
