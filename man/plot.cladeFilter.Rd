% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.cladeFilter.R
\name{plot.cladeFilter}
\alias{plot.cladeFilter}
\title{Plot \code{cladeFilter} Result in \code{trees} for One Grouping}
\usage{
\method{plot}{cladeFilter}(x, node.adj = c(0.5, 0.5), ...)
}
\arguments{
\item{x}{an object from "\code{cladeFilter}".}

\item{node.adj}{a numeric vector specifying the \code{adj} of node labels. Defaults
to \code{c(0.5, 0.5)}. Also see \code{plot.phylo}.}

\item{...}{other arguments to be passed to \code{plot} or to \code{plot.phylo}.}
}
\description{
Plot \code{trees} with tips colored for members of a grouping suggested to be
remained (green) or to be removed (red).
}
\seealso{
\code{\link{cladeFilter}}, \code{\link[ape]{plot.phylo}}
}
