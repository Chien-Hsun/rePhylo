
#' @title Summarize the Bootstrape Values of Backbone Nodes
#'
#' @description This function defines backbone nodes based on \code{taxa}
#' and \code{level}, summarizes the boostrape supports of these nodes,
#' and gives the order of \code{trees} by scoring the the distribution
#' of the collected bp values.
#'
#' @param trees a list of objects of class "\code{phylo}".
#'     A list of test trees such as single gene trees.
#' @param taxa a \code{data.frame} with tips and groupings. The first column
#'     must be the tip names in \code{trees}. Other columns specify the
#'     their groupings. Can have some tips left as no grouping.
#' @param level can be a character as one of the column names of \code{taxa},
#'     or a numeric specifying which column is used as grouping.
#'     Defaults as \code{NULL} and use the second column of \code{taxa}.
#' @param plot a logical specifying whether to generate violin plots
#'     for bootstrape values of backbone nodes of each of the \code{trees}.
#' @return the returned list contains 5 parts:
#'     (1) a \code{summary} table showing the quantiles and numbers
#'     of bp values >= 50 or 70 bp, and the score
#' @export
#' @details \code{backboneBP} first finds and uses only the trees meet all constraints
#'     (groupings) in \code{taxa} and \code{level} to identify backbone nodes. This function
#'     returns a list of four: (1) \code{summary}: main result, a table contains statistics
#'     and scores of \code{trees} that meets all constriants; (2) \code{data}: a list of modified
#'     \code{trees} that meets all constraints, to have bootstrape values on backbone nodes
#'     and \code{NA} on other nodes; (3) \code{plot.data}: a \code{data.frame} that can be as input
#'     to \code{ggplot2} to produce the violin plot; (4) \code{concor.trees}: statistics of
#'     the number of all \code{trees} meet constraints, including those meets partial or
#'     no constraints (groupings). This last result can present the consistencies of
#'     constraints in \code{trees}.
#' @import ggplot2
#' @importFrom ape root.phylo
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_violin
#' @importFrom ggplot2 geom_boxplot
#' @importFrom ggplot2 geom_dotplot
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 element_text
#' @importFrom stats quantile
#' @importFrom stats reorder
#' @examples
#' data(Brassidata)
#' trees <- Brassidata$trees
#' taxa <- Brassidata$taxaTable
#' ubb <- backboneBP(trees, taxa = taxa, level=2, plot = FALSE)
#' # returns NULL when there is no so-called "backbone nodes"
#' # based on the given groupings


# to note:
# concor.trees will perform in backboneBP , with its result shown in the object of backboneBP

backboneBP <- function(trees, taxa, level = NULL,
                        plot = TRUE){

  #################################################################
  # check inputs

  clas <- lapply(trees, function(z) inherits(z, "phylo"))
  clas <- unlist(clas)
  if(!all(clas))
    stop("Trees should be a list of objects of class \"phylo\".")

  treenames <- names(trees)
  if(is.null(treenames)){
    treenames <- c(1:length(trees))
  }

  # check taxa and level
  if(!inherits(taxa, "data.frame"))
    stop("Taxa should be an object of class \"data.frame\".")
  if(ncol(taxa) < 2)
    stop("Taxa should have at least two columns.")

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


  # performing concor.trees
  concor.tr <- concor.trees(trees = trees, taxa = taxa, level = level)

  const <- colnames(concor.tr$summary$constraints)
  score <- concor.tr[[1]][, "score"]
  # then to get the trees meets all clade constraints == score as 0
  # (if trees don't meet all constraints , then score will be a minus value)
  w <- which(score == 0); length(w)
  if(length(w) == 0){
    stop("No tree meets all groupings when performing concor.trees.")
  } else {
    if(length(w) != length(trees)){
      lenw <- length(w)
      cat(lenw, "trees fulfill all the groupings and are passed to the analysis.\n")
    }

    goodtrees <- vector("list", length(w))
    for(t in 1:length(w)){
      tt <- w[t]
      goodtrees[[t]] <- trees[[tt]]
    }
    treenames <- names(trees)
    names(goodtrees) <- treenames[w]
    #    trees <- goodtrees
  }
  # then the trees are those meets all the required constraints


  # 0724 to here
  ########################################################################################
  # TO DO:
  ########################################################################################
  # 0710
  # because the trees might not be rooted , so in some constraint that the member is as root
  # all the bp values is changed to NA by the original rationale
  # (first make a bp matrix , then change bp to NA for des nodes of each constraint)
  # MAYBE using the inverse idea , that to get the values as ancestors of the mrca of constraints
  # and then for those with all NA or the mrca == this is used as root
  # then root the tree again and just analyze this constraint again and include its values
  # THIS IDEA of positively grep the values seems applicable to all other functions!!!!????
  # and the way to know whether or not another rooting is required,
  # maybe is from whether the mrca of this constriant is "ROOT" in the test tree ???



    s.bpinfo <- goodtrees
  # to edit node.labels in trees directly
  # this copy step is required for repeating edit
  # for the nodel.label of the tree

    for(jj in 1:length(s.bpinfo)){
      s.tree <- s.bpinfo[[jj]]
      s.tree <- ape::root.phylo(phy = s.tree,
                                outgroup = s.tree$tip.label[1],
                                resolve.root = TRUE,
                                edgelabel = TRUE)
      s.tree<-rewrite.tree(s.tree)

      # doing level
      for(c in 1:length(const)){
        constr <- const[c]

        mrcalist <- .getmnode(constr = constr, tree = s.tree,
                             taxa = taxa, level = level,
                             bp = TRUE)
        # from is for bp values in "root.dist"

        if(!is.null(mrcalist)){
          s.tree <- select.backbone(tre = s.tree, mrcalist)
        }

      } # for c in const , end of doing level

      s.bpinfo[[jj]] <- s.tree

    } # for each tree (in s.bpinfo)
    names(s.bpinfo)<-names(goodtrees)


  s.bpinfo2 <- lapply(s.bpinfo, function(x) return(x$node.label[!is.na(x$node.label)]))
  # here the added "Root" and "" node.labels are set as NA , so is removed in the above line
  statl <- s.bpinfo2
  for(i in 1:length(s.bpinfo2)){
    sbp <- as.numeric(s.bpinfo2[[i]])
    ss <- stats::quantile(sbp,na.rm = TRUE)
    # quantile reflect the distribution of bp values
    # thus if higher values, indicating better supports for more nodes
    ss2 <- t(matrix(ss))
    colnames(ss2) <- names(ss)

    cc1 <- length(which(sbp >= 70))
    cc2 <- length(which(sbp <= 50))
    d <- data.frame(cc1, cc2, stringsAsFactors = FALSE)

    names(d) <- c("<=50bp", ">=70bp")
    dat <- cbind(ss2, d)

    statl[[i]] <- dat
  }


  # to make a stat table , with ordering for the goodness of trees
  cname <- colnames(statl[[1]])
  stab <- matrix(unlist(statl), nrow = length(statl), byrow = TRUE)
  stat <- data.frame(tree = names(s.bpinfo), stab)
  colnames(stat) <- c("tree", cname)
  #  head(stat);nrow(stat)

  # Description for the indicators for scoring:
  # for all indicators, the higher values are the better
  # for example:
  # "<=50bp", ">=70bp": indicate more node have bp values higher than 50 or 70 bp
  # "25%", "50%": indicate higher bp values in general
  # NOTE: in such case almost all indicators have some importance,
  # such as when 100% is not 100 means there is no bp 100 for any node!!

  # 2019.05.19
  # scoring by directly sum all the values in stat table
  # that is, the distribution of 0, 25, 50, 75, 100% of all backbone bp values
  # and also the number of nodes receiving bp >= 50 or 70
  score <- as.integer(apply(stat[,2:ncol(stat)], 1, sum))
  stat <- cbind(stat, score); stat
  stat <- stat[order(score, decreasing = TRUE), ]
  #  head(stat)

  sname <- names(s.bpinfo2)
  for(i in 1:nrow(stat)){
    tname <- stat[i, 1]
    w <- which(sname == tname)
    sbp <- as.numeric(s.bpinfo2[[w]])
    t <- rep(tname, length(sbp))
    ii <- rep(i, length(t))
    d <- data.frame(tree = t, or = ii, bp = sbp, stringsAsFactors = FALSE)
    if(i == 1)
      dd <- d
    else
      dd <- rbind(dd, d)
  }
  dd <- data.frame(dd)
  # there are some NAs by empty node.label in unroot trees
  dd <- dd[!is.na(dd[ ,"bp"]), ]
  # dd is the matrix for plot

  if(nrow(dd) == 0){

    g <- NULL

  } else {

    # to make violin plot with dots
    g <-
      ggplot2::ggplot(data = dd,
              ggplot2::aes_string(x = stats::reorder(x = dd[ ,"tree"],
                                  X = dd[ ,"or"]), y = "bp")) +
      ggplot2::geom_violin(ggplot2::aes_string(fill = "tree"),
                           trim = TRUE, show.legend = FALSE) +
      ggplot2::geom_dotplot(binaxis = 'y' , stackdir = 'center' ,
                            dotsize = .2) +
      #    ggplot2::geom_jitter(shape=16) +
      # , position=ggplot2::position_jitter(0.2)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                         hjust = 1)) +
      ggplot2::coord_flip() + ggplot2::geom_boxplot(width = 0.1) +
      ggplot2::xlab("Tree name")

    if(plot){
      print(g)
    }

  }

  resl <- list(summary = stat, data = s.bpinfo, plot.data = dd,
               ggplot = g, concor.trees = concor.tr)

  return(resl)
} # end of function





