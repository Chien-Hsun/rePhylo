#' @title Analyze Concordances on Trees with All Groupings
#'
#' @description An internal function of \code{backboneBP}. To summarize the concordances
#' with all groupings of each of \code{trees}, and returns trees that meets all groupings
#' to the main analysis of \code{backboneBP}.
#'
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of test trees such as single gene trees.
#' @param taxa a \code{data.frame} with tips and groupings. The first column must be
#'     the names of all tips in \code{trees}. Other columns specify the groupings.
#' @param level a character specifying the level used as grouping. Should be one of the
#'     column names of \code{taxa}. Defaults as \code{NULL} and use the second column
#'     of \code{taxa}.
#' @export
#' @seealso \code{\link{backboneBP}}
#' @examples
#' \dontrun{
#' data(Brassidata)
#' trees <- Brassidata$trees
#' taxa <- Brassidata$taxaTable
#' uconcort <- concor.trees(trees, taxa, level = 2)
#' }



# get the concordance results by trees, with specified nodes or clades (groupings)
# (for each of the trees, are they meet the requirement of all constraints?
# if not, how many meets all of them? and how many constraints are others meet?
# and, can take only those meet all requiremnts (can be optional).)

# "level" here is inherited from "backboneBP"
concor.trees <- function(trees, taxa, level = NULL){

  # TO DO:
  # 3. maybe also add a "bp" option???; can see internal.concor.node.R

  ########################################################################
  # check inputs
  if(!inherits(trees, "list"))
    stop("Input should be a list of objects of class \"phylo\".")

  t <- unlist(lapply(trees, function(x) inherits(x, "phylo")))
  if(!all(t))
    stop("Input should be a list of objects of class \"phylo\".")

  treenames <- names(trees)
  if(is.null(treenames)){
    treenames <- c(1:length(trees))
  }

  if(!inherits(taxa, "data.frame"))
    stop("Taxa should be an object of class \"data.frame\".")

  if(ncol(taxa) < 2)
    stop("Taxa should have at least two columns.")

  # "level" here is inherited from "backboneBP" and will be the ID in column number of taxa table
  # but considering the possible requirement for using concor.trees alone
  # still to check level
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

  #################################################################################
  # check if tip labels are all matched
  temp <- .match.tiplabels(trees = trees, taxa = taxa)
  trees <- temp$trees
  taxa <- temp$taxa

  ########################################################################
  # main content

  const <- taxa[ ,level]
  const <- const[!duplicated(const)]; const; length(const)

    checkconcor <- .check.concor(trees, taxa, level) # internal.backbp
    # as a result, if we want to get the trees with all the tips
    # grouped well as the constraints
    # then there is no tree can be used.
    # Before discuss with Hsun, I am going to write script
    # to just "remove" these separating tips
    # (all tips for that constraint) in the test trees



  # to get the number of constraints in each test trees
    # (in case there is missing constraint in some trees)
  allcon <- function(x, taxa, level){
    testtips <- x$tip.label
    alltips <- as.character(taxa[ ,1])

    w <- which(alltips %in% testtips)
    wc <- as.character(taxa[w, level])
    wc <- wc[!duplicated(wc)]

    wwc <- rep(FALSE, length(wc))
    for(i in 1:length(wc)){
      cc <- wc[i]
      tt <- as.character(taxa[taxa[ ,level] == cc,1]);tt
      tt <- tt[!is.na(match(tt, testtips))];tt
      if(length(tt) > 1){
        wwc[i] <- TRUE
      }
    }
    wc <- wc[wwc]

    return(wc)
  }
  allconst <- lapply(trees, function(x) allcon(x, taxa, level))
  # as denominators of constraints of each tree

  taball <- table((unlist(allconst)))
  # first to summarize the .concon result for each tree
  treescore <- unlist(lapply(checkconcor, length))
  # now treescore shows that each tree meets how many constraints
  # this could be one of the scoring
  ts <- treescore-taball
  sumtab1 <- rbind(treescore, taball, ts)
  row.names(sumtab1) <- c("concor", "all", "score")

  # to get the information of concor number of clades in each tree
  score <- rep(NA, length(trees))

  .temp <- function(x, i){
    ccc <- 0
    if(i %in% x)
      ccc <- 1
    return(ccc)
  }

  for(i in 1:length(trees)){
    ss <- lapply(checkconcor, function(x) .temp(x, i))
    score[i] <- sum(unlist(ss))
  }
  aa <- unlist(lapply(allconst, length))
  mm <- score-aa
  sumtab2 <- table(mm)
  sumtab2 <- t(matrix(as.numeric(sumtab2)))
  colnames(sumtab2) <- names(table(mm))
  row.names(sumtab2) <- "n.trees"
  sumlist <- list(constraints = sumtab1, concordances = sumtab2)

  dd <- data.frame(concor = score, allconst = aa, score = mm, stringsAsFactors = FALSE)


  ########################################################################

  result <- list(result = dd, summary = sumlist)
  # the new result contains a table with tree names and the number of concor for clades (constraints)
  # and a "sum" table of how many trees have groupings for each clade

  # in the result list, results are as the list of names of goodtrees and badtrees
  # score are the number of fulfilling the requirement of constraints for each tree
  # sum is the table of how many trees have groupings for each clade

  class(result) <- "concor.trees"
  return(result)
} # end of function concor.trees









################################################################################
# internal function
#' @title .check.concor
#' @description Internal function of \code{youshu}
#' @param tre tre
#' @param taxa taxa
#' @param level level
#' @export
.check.concor <- function(tre, taxa, level){
  # This is the original version, which applies for multiple trees with the same taxa and const information

  ##############################################################
  # some checks
  if(!inherits(taxa, "data.frame"))
    stop("Taxa table should be an object of class \"data.frame\".")
  t <- lapply(tre, function(x) inherits(x, "phylo"))
  if(!all(unlist(t)) | !inherits(tre, "list"))
    stop("Trees should be a list containing objects of class \"phylo\".")

  ##############################################################
  # main content
  const <- taxa[ ,level]
  const <- const[!duplicated(const)]; const; length(const)

  judgeTrees <- vector("list", length(const))
  names(judgeTrees) <- const

  tiplist <- as.list(const)
  names(tiplist) <- const

  for(c in 1:length(const)){
    constr <- const[c]
    tt <- taxa[ ,level] == constr
    w <- taxa[which(tt),1]
    tiplist[[c]] <- as.character(w)
  }

  tiplist2 <- tiplist
  ww <- c(1:length(tre))
  for(ll in 1:length(tiplist)){
    tiplist2[[ll]] <- list(tips = tiplist[[ll]], trees.w.node = ww)
  }


  # copy and edited from concor.node.unrooted, in order to avoid the use of ref tree
  m.concor.node <- function(tiplist2, trees, taxa){

    ##############################################################
#    testTree<-trees
#    testSet<-c(1:length(trees))
    labNA<-function(xx){
      xx$node.label<-as.character(as.numeric(xx$node.label))
      return(xx)
    }
    testTree<-lapply(trees,function(xx) labNA(xx))

    con_result <- lapply(tiplist2, function(x)
      .concon(x, testTree = testTree, setBP = 0, reftips = as.character(taxa[ ,1]),
             from = "concor.trees"))
    ##############################################################

    cr <- lapply(con_result, function(xx) .conRes(xx))

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

    # end of editing concor.node.unrooted
    #################################################################################################

    return(conTrees)
  } # end of m.concor.node function


  conTrees <- m.concor.node(tiplist2 = tiplist2, trees = tre, taxa = taxa)

  treenames <- names(tre)
  for(c in 1:length(conTrees)){
    goodt <- conTrees[[c]]
    wg <- treenames %in% goodt
    judgeTrees[[c]] <- which(wg)
  }

  return(judgeTrees)
} # end of ".check.concor"



