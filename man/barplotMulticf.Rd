% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barplotMulticf.R
\name{barplotMulticf}
\alias{barplotMulticf}
\title{Summarize \code{cladeFilter} Result in Barplot}
\usage{
barplotMulticf(
  x,
  type = c("groupings", "genetrees", "species"),
  reorder = "up",
  plot = TRUE,
  pdf = FALSE,
  pwd = c(7, 10, 10)
)
}
\arguments{
\item{x}{an object from "\code{summarizeMulticf}".}

\item{type}{c("\code{groupings}", "\code{genetrees}" or "\code{species}"),
specifying the type of the barplots. Defaults to all.}

\item{reorder}{one of "\code{up}", "\code{down}" or "\code{none}", specifying whether to reorder the barplots. Defaults to "\code{up}".}

\item{plot}{a logical specifying whether to plot barplots. Defaults to \code{TRUE}.}

\item{pdf}{a logical specifying whether to export plots into pdf files.
Defaults to \code{FALSE}.}

\item{pwd}{Numeric specifying the width of pdf files. Should be in the same length of "\code{type}".
Defaults to c(7,10,10).}
}
\description{
An internal function of \code{summarizeMulticf} to produce summary barplots.
}
\examples{
data(Brassidata)
trees <- Brassidata$trees
taxa <- Brassidata$taxaTable
# perform cladeFilter analysis
ures <- cladeFilter(trees = trees, taxa = taxa, level = 2)
# summarize cladeFilter results with multiple groupings
usum <- summarizeMulticf(ures, trees, plot.tree = FALSE, write.table = FALSE)
# returns NULL when the grouping covers the whole tree
if(!is.null(usum)){
# to make barplots
utest <- barplotMulticf(x = usum, reorder="up", plot = TRUE, pdf=FALSE)
}
}
\seealso{
\code{\link{summarizeMulticf}}
}
