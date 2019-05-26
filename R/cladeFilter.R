#' @title Filtering Sequences of Paralogy or With Other Misleading Signals Based on Given Groupings
#'
#' @description This function judges whether or not members of a grouping defined by \code{taxa}
#' and \code{level} grouped in \code{trees}, and suggest sequences to be removed when they are not.
#' Rationales for the selection of excluding sequences see \code{Details}.
#'
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of test trees such as single gene trees.
#' @param taxa a \code{data.frame} with tips and groupings. The first column must be
#'     the names of all tips in \code{trees}. Other columns specify the groupings.
#' @param level a character specifying the level used as grouping. Should be one of the
#'     column names of \code{taxa}. Defaults as \code{NULL} and use the second column
#'     of \code{taxa}.
#' @param write.table a logical specifying whether to write out text files
#'     containing the number and names of sequences suggested to be removed
#'     by \code{cladeFilter} for each grouping.
#'
#' @details To determine a scenario as the best solution of which sequences to be removed
#'     when members of a grouping is separated into more than one sub-clade, \code{cladeFilter}
#'     first calculate for each subclade the number of tips of a grouping to be remained or
#'     removed in the tree, and choose the one with the highest benefit as [number of
#'     remaining tips minus that of removing tips] (\code{max.benefit}). If this number is the
#'     same in more than one subclade, then \code{cladeFilter} chooses the one removing the
#'     least number of tips (\code{min.rmtips}). If both \code{max.benefit} and \code{min.rmtips}
#'     in more than one subclade are the same, \code{cladeFilter} searches the lengths of branches
#'     subtending to these subclades, remains the one with shortest branch length and
#'     removes the others. If the branch lengths are still the same for more than one
#'     subclade, \code{cladeFilter} then suggests to remove all tips of this grouping
#'     because the possibilities for these subclades to be paralogs are similar.
#'
#' @return A list of cladeFilter result for all test trees. For a test tree, the result constains:
#' result    result of cladeFilter
#' removingTips    tip(s) violate the grouping
#' remainingTips    tip(s) group
#' mrcaInTree    tips in the subclade of the MRCA node of the grouping in the test tree
#' node    MRCA node for the grouping
#' subtree    the extracted subclade of the grouping
#'
#' @export
#' @importFrom ape extract.clade
#' @importFrom ape getMRCA
#' @importFrom ape is.rooted
#' @importFrom phangorn Descendants
#' @importFrom phangorn Ancestors
#' @importFrom utils write.table
#' @seealso \code{\link{summarizeMulticf}}, \code{\link{barplotMulticf}}, \code{\link{plot.cladeFilter}}
#' @examples
#' \dontrun{
#' data(Brassidata)
#' trees <- Brassidata$trees
#' taxa <- Brassidata$taxaTable
#' # perform cladeFilter analysis
#' ures <- cladeFilter(trees = trees, taxa = taxa, level = 2, write.table = FALSE)
#' # summarize cladeFilter results with multiple groupings
#' usum <- summarizeMulticf(ures, trees, plot.tree = FALSE, write.table = FALSE)
#' # to make barplots
#' utest <- barplotMulticf(x = usum, reorder="up", plot = TRUE, pdf=FALSE)
#'
#' }



cladeFilter <- function(trees, taxa, level = NULL, write.table = TRUE){
  # removing "warn = FALSE"

  warn <- TRUE
  #  library(ape) # required
  #  library(phangorn) # required

  clas <- unlist(lapply(trees, function(z) inherits(z, "phylo")))
  if(!inherits(trees, "list") | !all(clas))
    stop("Trees should be a list of objects of class \"phylo\".")

  fns <- names(trees)
  if(is.null(fns)){
    fns <- c(1:length(trees))
  }

  if(!inherits(taxa, "data.frame"))
    stop("Taxa should be an object of class \"data.frame\".")

  # to check level
  if(is.null(level)){
    # if there is no specified level for taxa table,
    # will use the second column directly
    level<-2
  } else {
    # if there is level, considering them (1) as number (2) as character
    if(inherits(level, "character")){
      levels<-colnames(taxa)
      level <- match.arg(level, choices = levels, several.ok = FALSE)
      level <- which(colnames(taxa) == level)
    }
    # if level is numeric, don't do checks
  }


  ########################################################################
  # check if tip labels are all matched
  temp <- .match.tiplabels(trees = trees, taxa = taxa)
  trees <- temp$trees
  taxa <- temp$taxa


  ########################################################################
  # to make g list from taxa and level
  const <- taxa[ ,level]
  const <- const[!duplicated(const)]
  g <- vector("list", length(const))
  for(i in 1:length(const)){
    name <- const[i]
    grouping <- taxa[taxa[ ,2] == name,1]
    grouping <- data.frame(grouping, stringsAsFactors = FALSE)
    r <- list(grouping = grouping, name = name)
    g[[i]] <- r
  }
  names(g) <- const


  ###############################################################
  # main codes are using ".cladefilter"
  # will root the trees within ".cladeFilter" based on group if trees are not all rooted
  res <- lapply(g, function(xx) .cladefilter(trees, xx, warn, fns))
  names(res) <- const


  if(write.table){

    for(i in 1:length(res)){
      gdat <- res[[i]]
      treenames <- names(gdat)
      if(is.null(const)) gname <- NULL
      else gname <- const[i]

      removeNum <- removeTable <- NULL
      for(j in 1:length(gdat)){
        tdat <- gdat[[j]]
        removeTips <- tdat[[2]]
        if(!is.null(removeTips)){
          len <- length(removeTips)
          if(len > 0){
            if(is.null(removeNum)){
              removeTable <- cbind(treefile = treenames[j], removeTips = removeTips)
              removeNum <- cbind(treefile = treenames[j], removeNum = len)
            } else {
              temp <- cbind(treefile = treenames[j], removeTips = removeTips)
              removeTable <- rbind(removeTable, temp)
              temp <- cbind(treefile = treenames[j], removeNum = len)
              removeNum <- rbind(removeNum, temp)
            }
          }
        }
      }

      if(is.null(gname)){
        if(!is.null(removeTable)){
          fname <- "remove_Table.out"
          utils::write.table(removeTable, file = fname, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
          fname <- "remove_number.out"
          utils::write.table(removeNum, file = fname, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
      } else {
        if(!is.null(removeTable)){
          fname <- paste("remove_Table_", gname, ".out", sep = "")
          utils::write.table(removeTable, file = fname, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
          fname <- paste("remove_number_", gname, ".out", sep = "")
          utils::write.table(removeNum, file = fname, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
        }
      }
    }
  }

  return(res)

} # end of function



#############################################################################
#' @title .cladeFilter
#' @description Internal function of \code{BetterTree}
#' @param trees trees
#' @param xx xx
#' @param warn warn
#' @param fns fns
#' @importFrom ape is.rooted
#' @export
.cladefilter <- function(trees, xx, warn, fns){

  fres <- vector("list", length(trees))
  group <- xx$grouping
  err <- xx$name
  isroot <- unlist(lapply(trees, ape::is.rooted))
  if(!all(isroot)){
    trees <- lapply(trees,function(zz) root.dist(z = zz, tar = group[,1]))
  }

  for(z in 1:length(trees)){

    tree <- trees[[z]]
    fname <- fns[z]

    # to format the mode of constraint
    if(!inherits(group, "data.frame"))
      stop("Grouping should be an object of class \"data.frame\".")

    # consider if there is only 1 tip in the constraint, there is no need to perform filtering
    # this is possible when using manually defined constraints
    if(nrow(group) == 1)
      stop("Only one tip in the grouping")

    if(is.null(tree)){
      # if there is no any tip of this group in the test tree
      # a NULL will be returned to trees by "root.dist"

      subtree <- NULL
      mrcaInTree <- NULL
      fnode <- NULL
      removeTips <- NULL
      untar_list <- NULL

    } else {

      # get only tips in tree of the MRCA file
      treeTips <- tree$tip.label
      matchTips <- match(group[ ,1], treeTips)
      mrcaInTree <- group[!is.na(matchTips),1]; length(mrcaInTree)
      mrcaInTree <- as.character(mrcaInTree); length(mrcaInTree)

      # if there is no tip or only one tip in the gene tree belongs to the constraint
      # then there is no sequence to be removed
      # Although in the given constraint there are always tips within
      # but it is uncertain in the individual gene trees
      if(length(mrcaInTree) == 0 | length(mrcaInTree) == 1){

        subtree <- NULL
        mrcaInTree <- NULL
        fnode <- NULL
        removeTips <- NULL
        untar_list <- NULL

      } # for length(mrcaInTree) == 0 or 1

      #####################################################################
      # 0524 to HERE
      # all the codes are checked, only left the result object to be sure
      # for below (before 214) (above are checked)
      # PLEASE see important notes on BLACK note book
      #####################################################################
      # newly added this option "judgeType", for two kinds of ways to decide which tip to be remvoved
      # 1. by "branch length", that even there are only 2 tips, if they are not grouped together,
      # compare their branch lengths and remove the one with longer branch
      # this is also applicable to subclades with multiple tips
      # 2. by "tip number", that if there are 2 tips and NOT grouped together,
      # then the two tips are both removed,
      # BUT for multiple tips for subclades, the "benefit" is based to select which subclade to be removed
      #####################################################################
      # if there is two tips
      if(length(mrcaInTree) == 2){
        # to see if they are together or not
        AncNode <- ape::getMRCA(tree, mrcaInTree)
        subtree <- ape::extract.clade(tree, AncNode)
        subtips <- subtree$tip.label
        matchSub <- match(subtips, mrcaInTree)
        jumpInTips <- subtips[is.na(matchSub)]
        # if the two tips grouped together
        if(length(jumpInTips) == 0){
          sub_root <- ape::getMRCA(subtree, mrcaInTree)
          fnode <- sub_root
          removeTips <- NULL
          untar_list <- NULL
        } else {
          # if the two tips separates from each other
          # remove the one with long branch
          bls <- .compare.edge.len(tree, subtips)
          ww <- which(bls == max(bls))
          if(length(ww) == 1){
            removeTips <- subtips[ww]
            fnode <- ape::getMRCA(subtree, mrcaInTree)
            untar_list <- list(num_wanted = 1, removeTips = removeTips, num_of_rm = 1, num_of_benefit = 0)
            untar_list <- list(untar_list) # !!!!!!!!!!!!!!!!!
            names(untar_list) <- fnode # !!!!!!!!!!!!!!!!!!!!!
          } else {
            # in case the lengths are the same
            # such case should be rare for ML trees
            # but will be common for coalescence trees !!
            if(warn)
              cat(z, ": ", fname, " in ", err, " grouping will have 2 tips separated with the same branch lengths. Both tips will be removed.\n", sep="")
            removeTips <- mrcaInTree
            fnode <- ape::getMRCA(subtree, mrcaInTree)
            untar_list <- list(num_wanted = 0, removeTips = removeTips, num_of_rm = 2, num_of_benefit=-2)
            untar_list <- list(untar_list) # !!!!!!!!!!!!!!!!!
            names(untar_list) <- fnode # !!!!!!!!!!!!!!!!!!!!!
          }
        }
      } # for length(mrcaInTree) == 2

      # required result objects:
      # (all the same=> fns, mrcaInTree), (most the same=> subtree),
      # each different=> removeTips, fnode, untar_list
      if(length(mrcaInTree) > 2){
        AncNode <- ape::getMRCA(tree, mrcaInTree)
        subtree <- ape::extract.clade(tree, AncNode)
        subtips <- subtree$tip.label
        matchSub <- match(subtips, mrcaInTree)
        jumpInTips <- subtips[is.na(matchSub)]; length(jumpInTips)

        if(length(jumpInTips) == 0){ # all mrca tips well grouped together
          # need to check the required objects for here
          removeTips <- NULL
          untar_list <- NULL

          sub_root <- ape::getMRCA(subtree, mrcaInTree)
          fnode <- sub_root
        }

        if(length(jumpInTips) == 1){ # if there is only 1 tip inserted, just remove it
          sub_root <- ape::getMRCA(subtree, mrcaInTree)
          fnode <- sub_root

          removeTips <- jumpInTips
          wantedTips <- subtips[!is.na(matchSub)]
          lr <- length(removeTips)
          lw <- length(wantedTips)
          benefit <- lw-lr
          untar_list <- list(lw, removeTips, lr, benefit)
          names(untar_list) <- c("num_wanted", "removeTips", "num_of_rm", "num_of_benefit")
          untar_list <- list(untar_list) # !!!!!!!!!!!!!!!!!
          names(untar_list) <- fnode # !!!!!!!!!!!!!!!!!!!!!
        }

        if(length(jumpInTips) > 1){
          # if there are several jumpInTips, the no. of jump tips might be larger than one of the subclade,
          # so it is possible to remove one of the subclades rather than the inserted tips (if inserted are so many)
          # we need to get information of all subclades and the inserted groups
          # and compare to know which case removes the least no. of tips

          # to find the split clades
          sub_root <- ape::getMRCA(subtree, mrcaInTree)
          allDes <- phangorn::Descendants(subtree, sub_root, type = "all")
          allDes <- allDes[allDes > length(subtree$tip.label)] # remain only nodes but not tips
          subclade_id <- 0
          for(i in 1:length(allDes)){
            node <- allDes[i]

            split_clade <- ape::extract.clade(subtree, node)
            split_tips <- split_clade$tip.label; length(split_tips)
            split_match <- match(split_tips, mrcaInTree)
            na <- is.na(split_match)
            if(!any(na)){
              subclade_id <- c(subclade_id, node)
            }
          }
          subclade_id <- subclade_id[-1]
          # here are all nodes that the clade contains only wanted species (will have nested nodes)

          if(length(subclade_id) == 0){
            # this is possible since the root is not included in allDes
            # and for any node there should be at least 2 tips
            # if length == 0, means all the mrca tips are separated without group into any clade
            # then just removing all these mrca tips
            sub_root <- ape::getMRCA(subtree, mrcaInTree)
            fnode <- sub_root

            removeTips <- mrcaInTree
            lr <- length(removeTips)
            lw <- 0
            benefit <- lw-lr
            untar_list <- list(lw, removeTips, lr, benefit)
            names(untar_list) <- c("num_wanted", "removeTips", "num_of_rm", "num_of_benefit")
            untar_list <- list(untar_list) # !!!!!!!!!!!!!!!!!
            names(untar_list) <- fnode # !!!!!!!!!!!!!!!!!!!!!
          } else {

            # to get the number of node of all the subclades
            # removing the nested nodes and remaining only the deepest MRCA of each subclade
            targetNode <- 0
            for(i in 1:length(subclade_id)){
              node <- subclade_id[i]; node
              anc <- phangorn::Ancestors(subtree, node, type = "parent"); anc
              if(!anc %in% subclade_id){
                targetNode <- c(targetNode, node)
              }
            }
            targetNode <- targetNode[-1]
            # this is the number of the subclades and the mrca ID
            # if there is only one targetNode, the process should still be the same as multiple targetNodes
            # that to expand the targetNode and to see which is with the best benefit
            # it is just to keep the only one node in tar_list2 (for example)

            # if more than 2 target nodes, then see if any subclades can be combined to a larger clade with a few tips jumped in
            # to know which of them can combined as one clade with insertion of few (that we can accepted) species
            tar_list <- vector("list", length(targetNode))
            node <- 0 # just to reset
            for(i in 1:length(targetNode)){
              node <- targetNode[i]
              anc <- phangorn::Ancestors(subtree, node, type = "all"); anc

              temp <- c(0, 0, 0)
              node2 <- 0
              for(j in 1:length(anc)){
                node2 <- anc[j]
                sclade <- ape::extract.clade(subtree, node2)
                stips <- sclade$tip.label; stips
                match <- match(stips, mrcaInTree); match
                no.leng <- length(stips[is.na(match)]); no.leng
                # the no. of added "wrong tips" in this anc node
                yes.leng <- length(stips[!is.na(match)]); yes.leng
                # the no. of added "wanted tips" in this anc node

                temp <- rbind(temp, c(node2, no.leng, yes.leng))
              }
              temp <- temp[-1, ]
              # to make temp with one row as matrix
              # for following codes easier
              if(inherits(temp, "numeric")){
                temp <- t(as.matrix(temp))
              }
              colnames(temp) <- c("node", "removed", "wanted")

              tar_list[[i]] <- temp
            }
            names(tar_list) <- targetNode#; tar_list; str(tar_list)
            # tar_list contains information for each ancestral node assume to be combined as a larger clade
            # and the no. of added removed/wanted tips to see the benefit of this combination


            # tar_mat is to judge whether the combination is good
            tar_mat <- vector("list", length(tar_list))
            for(i in 1:length(tar_list)){
              tnode <- targetNode[i]
              sclade <- ape::extract.clade(subtree, tnode)
              stips <- sclade$tip.label
              match <- match(stips, mrcaInTree)
              no.leng <- length(stips[is.na(match)]); no.leng # "wrong tips" for the original targetNode
              yes.leng <- length(stips[!is.na(match)]); yes.leng # "wanted tips" for the original targetNode
              t.benefit <- (yes.leng)-(no.leng); t.benefit

              t <- tar_list[[i]]
              node <- 0
              rm <- tar_list[[i]][ ,2]
              add <- tar_list[[i]][ ,3]

              benefit <- add-rm; benefit
              benefit <- benefit-(t.benefit); benefit
              # only want the combination which will cause an increase in the (wanted-removed) tips

              if(!any(benefit > 0)){
                node <- 0; node
                combine <- FALSE
              } else {
                which <- which(benefit == max(benefit)); which
                which <- which[length(which)]; which # use the deepest node
                if(which > 1){
                  node <- t[max(which),1]; node
                  combine <- TRUE
                } else {
                  node <- 0; node
                  combine <- FALSE
                }
              }
              temp <- list(node, combine)
              tar_mat[[i]] <- temp
            }
            # tar_mat should now be the best subclade solution to be set for each targetNode

            tar_list2 <- tar_list
            for(i in 1:length(tar_list)){
              if(tar_mat[[i]][[2]]){
                tar_list2[[i]] <- tar_mat[[i]][[1]]
              } else {
                tar_list2[[i]] <- names(tar_list)[i]
              }
            }; tar_list2
            untar <- as.numeric(unlist(tar_list2)); untar
            untar <- untar[!duplicated(untar)]; untar
            # untar is now the mrca of final subclades (after combination/or not)


            # then to find wanted/removed tips of each subclade
            # to see which have highest "percentage"
            untar_list <- vector("list", length(untar))
            n.wanted <- numeric(length(untar))
            n.rm <- numeric(length(untar))
            n.benefit <- numeric(length(untar))
            node <- 0
            for(i in 1:length(untar)){
              node <- untar[i]

              sclade <- ape::extract.clade(subtree, node)
              stips <- sclade$tip.label # tips in that subclade
              match <- match(stips, mrcaInTree)
              rm_out_tips <- stips[is.na(match)]; rm_out_tips # who jumped in
              target_tips <- stips[!is.na(match)]; target_tips # remains in this clade
              n.wanted[i] <- length(target_tips)

              match <- match(mrcaInTree, stips)
              rm_in_tips <- mrcaInTree[is.na(match)]; rm_in_tips # mrca tips not in this subclade

              removeTips <- c(rm_in_tips, rm_out_tips) # all the removed tips
              n.rm[i] <- length(removeTips)
              n.benefit[i] <- n.wanted[i]-n.rm[i]

              untar_list[[i]] <- list(n.wanted[i], removeTips, n.rm[i], n.benefit[i])
              names(untar_list[[i]]) <- c("num_wanted", "removeTips", "num_of_rm", "num_of_benefit")
            }
            names(untar_list) <- untar

            # to find the subclade with highest percent of wanted tips when the clade is to be remained
            node <- 0
            which <- which(n.benefit == max(n.benefit)); which

            if(length(which) > 1){
              which <- which(n.rm == min(n.rm)); which
            }
            # if the n.benefit are the same in subclades, then use the one removing least tips

            # if there are more than one nodes having the same max n.benefit
            # then should choose the one with shortest branch length in its mrca node
            # but the possibility should be low for such case
            # so let the warning goes first and edit the code latter
            if(length(which) > 1){  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              if(warn)
                cat(paste("No. ", z, ". ", fname, " in ", err, " grouping will have more than one node with the same max_n.benefit and min_n.rm.\n", sep=""))

              # to judge the branch lengths of MRCA branch, and then choose the shortest one
              # first get the mrca node of these subclades
              allbr <- numeric(length(which))
              for(nn in 1:length(which)){
                ww <- which[nn]
                problemN <- targetNode[ww]
                brlen <- which(subtree$edge[ ,2] == problemN)
                allbr[nn] <- subtree$edge.length[brlen]
              }
              br_min <- which(allbr == min(allbr))
              if(length(br_min) == 1){
                which <- which[br_min]
              } else {
                if(warn)
                  cat("No. ", z, ". ", fname, " in ", err, " grouping even have the same branch lengths. All the subclades will be removed.\n",  sep="")

                removeTips <- mrcaInTree
                fnode <- sub_root
              }
            }

            if(!is.null(which)){ # when "which == NULL" indicates it goes to multiple br_max!!!!
              # then remove all the mrcaInTree because which the paralog is is hard to determine
              # and the sequences seems confusing
              fnode <- untar[which]; fnode

              sclade <- ape::extract.clade(subtree, fnode)
              stips <- sclade$tip.label
              match <- match(stips, mrcaInTree)
              rm_out_tips <- stips[is.na(match)]#; rm_out_tips

              match <- match(mrcaInTree, stips)
              rm_in_tips <- mrcaInTree[is.na(match)]#; rm_in_tips

              removeTips <- c(rm_in_tips, rm_out_tips)

            }

          }
        }
      }

      if(!is.null(removeTips)){
        # green is the remaining mrca tips and red is the removing tips
        green <- which(subtree$tip.label %in% mrcaInTree)
        red <- which(subtree$tip.label %in% removeTips)
        # must first do green then do red
        # becasue some mrca tips might be removing tips
        which <- which(green %in% red)
        if(length(which) != 0){
          wanted <- green[-which] # green[which] is the mrca tips in the separated/removing subclade
          wanted <- subtree$tip.label[wanted]
        } else {
          wanted <- subtree$tip.label[green]
        }
      } else {
        wanted <- mrcaInTree
      }

    } # end of if tree is not NULL (that there are tips for this group)

    fres[[z]] <- list(result = untar_list, removingTips = removeTips, remainingTips = wanted,
                      mrcaInTree = mrcaInTree, node = fnode, subtree = subtree)

  } # for loop of gene trees

  names(fres) <- fns
  class(fres) <- "cladeFilter"
  return(fres)
}

