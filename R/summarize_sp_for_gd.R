#' Summarize number of species supporting gd
#'
#' Summarize the number of gd supported by each species, 
#'     and the number of gd supported by how much species.
#'
#' @param phyto a data.frame of gd table of a specific phyto node.
#' @param pdf a logical specifying whether to export plots into pdf files.
#'     Defaults to \code{TRUE}.
#' @param pdfwid a numeric for width of pdf file. Defaults to 5.
#' @param pdflen optional. Length of pdf file. Defaults to \code{NULL}.
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
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @importFrom stats reorder

summarize_sp_for_GD <- 
  function(phyto, pdf = TRUE, pdfwid = 5, pdflen = NULL){
  
  pname <- phyto[1, 3]
  
  allgd <- unique(phyto[ , 2])
  res <- lapply(allgd, function(xx) 
    .sum_sp_gd(xx, phyto))
  
  # make res into table
  dat <- c()
  for(i in 1:length(allgd)){
    rr <- res[[i]]
    dat <- rbind(dat, rr)
  }
  dat <- data.frame(dat, stringsAsFactors = FALSE, 
                    row.names = NULL)
  colnames(dat) <- c("gd", "sp_num", "sp")
  # write out table
  wname <- paste("Summarize_sp_supporting_", 
                 pname, ".txt", sep="")
  write.table(dat, wname, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  
  # make plot-1: x == each sp, y == count of gd
  allsp <- dat[ , "sp"]
  allsp <- unlist(strsplit(allsp, split = ",", 
                           fixed = TRUE))
  sp <- allsp
  data <- data.frame(sp)
  
  # make plot
  gglist <- vector("list", 2)
  g1 <- ggplot2::ggplot(data, ggplot2::aes_string(x = "sp"))
  
  gg1 <- g1 + ggplot2::geom_bar(ggplot2::aes_string(fill = "sp"),
                                position = ggplot2::position_stack(reverse = TRUE),
                                show.legend = FALSE) +
    ggplot2::coord_flip() + 
    ggplot2::guides(fill=ggplot2::guide_legend(title = "")) + 
    # no title for color legend
    ggplot2::xlab(label = "sp")
  
  gglist[[1]] <- gg1
  
  if(pdf){
    if(is.null(pdflen)){
      # set pdf size
      pdflen <- unique(data[ , "sp"])
      pdflen <- length(pdflen)
      pdflen <- pdflen * 0.2
      if(pdflen < 3)
        pdflen <- 3
    }

    pdfname1 <- paste("Distribution_of_species_", 
                      pname, ".pdf",sep="")
    grDevices::pdf(pdfname1, pdfwid, pdflen)
    print(gg1)
    grDevices::dev.off()
  }
  
  # make plot-2: x == sp number, y == count of gd
  sp_num <- dat[ , "sp_num"]
  data <- data.frame(sp_num)
  
  # make plot
  g1 <- 
    ggplot2::ggplot(
      data, ggplot2::aes_string(
        x = stats::reorder(x = data[ , "sp_num"], 
                         X = rep(1, nrow(data)), 
                         FUN = function(x) sum(x))))
  gg1 <- g1 + ggplot2::geom_bar(
    ggplot2::aes_string(fill = "sp_num"), 
    position = ggplot2::position_stack(reverse = TRUE),
                                show.legend = FALSE) +
    ggplot2::coord_flip() + 
    ggplot2::guides(fill = ggplot2::guide_legend(title = "")) + 
    # no title for color legend
    ggplot2::xlab(label = "sp_num")
  
  gglist[[2]] <- gg1
  
  if(pdf){
    if(is.null(pdflen)){
      # set pdf size
      pdflen <- unique(data[ , "sp_num"])
      pdflen <- length(pdflen)
      pdflen <- pdflen * 0.2
      if(pdflen < 3)
        pdflen <- 3
    }
    
    pdfname2 <- paste("Distribution_species_num_", 
                      pname, ".pdf", sep="")
    grDevices::pdf(pdfname2, pdfwid, pdflen)
    print(gg1)
    grDevices::dev.off()
  }
  
  return(gglist)
}




# internal small function
#' Internal function
#'
#' Internal function of summarize_sp_for_GD
#'
#' @param xx a gd id
#' @param phyto a data.frame of gd table of a specific phyto node.
#' @export
#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar

.sum_sp_gd <- function(xx, phyto){
  gd <- xx
  w <- which(phyto[ , 2] == gd)
  sp <- phyto[w, 4]
  sp <- unique(sp)
  len <- length(sp)
  sp <- paste(sp, collapse = ",")
  rr <- c(gd, len, sp)
  return(rr)
}


