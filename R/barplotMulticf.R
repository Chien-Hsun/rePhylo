#' Summarize \code{cladeFilter} Result in Barplot
#'
#' An internal function of \code{summarizeMulticf} to produce summary barplots.
#'
#' @param x an object from "\code{summarizeMulticf}".
#' @param type c("\code{groupings}", "\code{genetrees}" or "\code{species}"),
#'     specifying the type of the barplots. Defaults to all.
#' @param reorder one of "\code{up}", "\code{down}" or "\code{none}", specifying whether to reorder the barplots. Defaults to "\code{up}".
#' @param plot a logical specifying whether to plot barplots. Defaults to \code{TRUE}.
#' @param pdf a logical specifying whether to export plots into pdf files.
#'     Defaults to \code{FALSE}.
#' @param pwd Numeric specifying the width of pdf files. Should be in the same length of "\code{type}".
#'     Defaults to c(7,10,10).
#' @export
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 position_stack
#' @importFrom ggplot2 guide_legend
#' @importFrom grDevices dev.off
#' @importFrom stats reorder
#' @seealso \code{\link{summarizeMulticf}}
#' @examples
#' data(Brassidata)
#' trees <- Brassidata$trees
#' taxa <- Brassidata$taxaTable
#' # perform cladeFilter analysis
#' ures <- cladeFilter(trees = trees, taxa = taxa, level = 2)
#' # summarize cladeFilter results with multiple groupings
#' usum <- summarizeMulticf(ures, trees, plot.tree = FALSE, write.table = FALSE)
#' # returns NULL when the grouping covers the whole tree
#' if(!is.null(usum)){
#' # to make barplots
#' utest <- barplotMulticf(x = usum, reorder="up", plot = TRUE, pdf=FALSE)
#' }


# 0801:
# include this function into "summarize.multiccf
# but still to export the function (just as suggested by "Writing R Extensions")
# 2019.05:
# pull out this function while re-using this package again
# in this way it will be more flexible

# to plot summarized barplots
# the result is a list of objects of class "ggplot"
# users can specify their preferences for plots by "+"
# if users want to plot from the result object (is a clear version), start from +geom_bar()

barplotMulticf <- function(x, type=c("groupings", "genetrees", "species"),
                    reorder="up", plot = TRUE, pdf = FALSE,
                    pwd=c(7,10,10)){

  if(!inherits(x, "data.frame"))
    stop("x should be an object of \"data.frame\".")

  if(ncol(x) != 3)
    stop("x should be a data frame with 3 columns.")

#  if(!missingArg(type)){
    len.t<-length(type)
    len.s<-length(pwd)
    if(len.t != len.s){
      stop("Arg \"pwd\" should be in the same length as \"type\"")
    }
#  }

  reorder<-match.arg(reorder,choices = c("up","down","none"),several.ok = FALSE)
  type<-match.arg(type,choices = c("groupings", "genetrees", "species"),several.ok = TRUE)

  cname <- colnames(x)
#  acname <- c("groupings", "genetrees", "species")
#  match <- match(cname, acname)
#  if(match[1] != 1 | match[2] != 2 | match[3] != 3)
#    warning("Please make sure the columns are in the order of \"groupings\", \"genetrees\", and \"species\".")

  # removed numbers, in genetrees, species (such as transcriptomes), and groupings
  data <- data.frame(x, row.names = c(1:nrow(x)))

  if(ncol(data) == 1){
    data <- as.data.frame(t(data));data
  }
  colnames(data) <- cname

  gplots<-as.list(c(1:length(type)))
  for(xx in 1:length(type)){
    ty<-type[xx]

    if(reorder == "up"){
      g1 <- ggplot2::ggplot(data, ggplot2::aes_string(x = stats::reorder(x=data[,ty], X=rep(1,nrow(data)), FUN=function(x) sum(x))))
    }
    if(reorder == "down"){
      g1 <- ggplot2::ggplot(data, ggplot2::aes_string(x = stats::reorder(x=data[,ty], X=rep(1,nrow(data)), FUN=function(x) -sum(x))))
    }
    if(reorder == "none"){
      g1 <- ggplot2::ggplot(data, ggplot2::aes_string(x = ty))
    }

    gg1 <- g1 + ggplot2::geom_bar(ggplot2::aes_string(fill = "groupings"),
                    position = ggplot2::position_stack(reverse = TRUE)) +
      ggplot2::coord_flip() + ggplot2::theme(legend.position = "top",
                           legend.box.just="right") +
      ggplot2::guides(fill=ggplot2::guide_legend(title="")) + # no title for color legend
      ggplot2::xlab(label = ty)

    gplots[[xx]]<-gg1
  } # for types
  names(gplots)<-type

  # return types
  if(plot & !pdf){
    print(gplots)
  }

  if(plot & pdf){
    for(xx in 1:length(type)){
      ty<-type[xx]

      # set pdf size
      len <- data[ ,ty]; length(len)
      len <- len[!duplicated(len)]
      len <- length(len)
      len <- len*0.2
      if(len < 4.5)
        len <- 4.5

      psize <- c(pwd[xx], len) # c(宽, 长)
#      g<-grep(x = ty,pattern = "groupings")
#      if(length(g) > 0){
#        psize <- c(pwd-2, len) # c(宽, 长)
#      }

      pdfname<-paste("barplot_",ty,".pdf",sep="")

      pdf(pdfname, psize[1], psize[2])
      print(gplots[[xx]])
      grDevices::dev.off()
    } # end of types in plot & pdf
  } # end of plot & pdf

  return(gplots)
}
