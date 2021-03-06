% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_GD_by_coverage.R
\name{filter_GD_by_coverage}
\alias{filter_GD_by_coverage}
\title{filters the GD counts by the species coverage of
    each of the duplicated subclades}
\usage{
filter_GD_by_coverage(
  phyto,
  node,
  trees,
  gdtable,
  sptree,
  split,
  tree_id_tab,
  sub_coverage = NULL,
  max.tip = NULL,
  up_one_node = FALSE,
  plot = FALSE,
  pdfwid = 10,
  pdflen = 30,
  mc.cores = 4
)
}
\arguments{
\item{phyto}{the phyto id for the analysis.}

\item{node}{a numeric indicating the target node id as in \code{ape}.
An alternative arg of \code{phyto}.}

\item{trees}{a list of objects of class "\code{phylo}".
A list of gene trees to be examined. Can be unrooted.}

\item{gdtable}{a \code{data.frame} of the gd table file provided by Tree2GD.}

\item{sptree}{the species tree used to do tree reconciliation.}

\item{split}{a character. The symbol used to separate species name and the 
sequence number in the gene family trees.}

\item{tree_id_tab}{a table for the tree id during Tree2GD and the name of the 
corresponding tree file. The first column is the id, and the second column
is the tree names.}

\item{sub_coverage}{a numeric specify the species coverage level. 
If <= 1, treated as percentage. If > 1, treated as number of species.}

\item{max.tip}{an integer specifying the max number of tips
for the coverage filtering. Applied only when sub_coverage 
is set as percent (<= 1). Used as a cutoff for coverage as 
N tips as for deep nodes having large clades.}

\item{up_one_node}{logical. Whether to choose tips of one node deeper as the 
closest outgroup, as a second step of investigation. Defaults to \code{FALSE}.}

\item{plot}{whether to plot subclades. Defaults to \code{FALSE}.}

\item{pdfwid}{the width of pdf file for subtrees. Defaults to \code{10}.}

\item{pdflen}{the length of pdf file for subtrees. Defaults to \code{30}.}

\item{mc.cores}{the number of cores used for mclapply. Defaults to \code{4}.}
}
\description{
\code{filter_GD_by_coverage} is to identify for each gd whether the two 
    duplicated subclades contains species number more than a given level.
}
\details{
returns a pdf plotting the gd clades, with gd pairs 
    in red and other species in the same clade in sptree in green.
    A file with prefix "Include_or_not_" indicates for each tree 
    and gd the requirement is fulfilled or not. A file with prefix
    "Trees_have_basal_tips_" shows the page number the tree fulfilled
    the requirement in the pdf file.
}
