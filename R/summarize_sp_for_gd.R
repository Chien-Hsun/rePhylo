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
#' @param plot.type characters to specify the type of plot. One of 
#'     "two", "sp", "sp_num". Can be more than one.
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
  function(phyto, pdf = TRUE, pdfwid = NULL, pdflen = NULL,
           plot.type = "two"){
  
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
  
  ress<-list()
  if(plot.type == "two"){
    #######################################################
    ## two-dimension plot
    #######################################################
    # prepare data frame for x-y dotplot
    .melt.sp.GD<-function(xx,dat){
      num<-dat[xx,2]
      sp<-dat[xx,3]
      sp<-unlist(strsplit(sp,split = ",",fixed = TRUE))
      num<-rep(num,length(sp))
      dd<-paste(num,sp,sep=" ")
      return(dd)
    }
    
    meltd<-lapply(c(1:nrow(dat)),function(xx) .melt.sp.GD(xx,dat))
    mdata<-c()
    for(i in 1:length(meltd)){
      dd<-meltd[[i]]
      mdata<-c(mdata,dd)
    }
    mtab<-table(mdata)
    size<-as.numeric(mtab)
    tmp<-names(mtab)
    tmp<-lapply(tmp,function(xx) 
      unlist(strsplit(xx,split = " ",fixed = TRUE)))
    sp_num<-unlist(lapply(tmp,function(xx) xx[1]))
    sp<-unlist(lapply(tmp,function(xx) xx[2]))
    mdata<-data.frame(sp_num,sp,size,stringsAsFactors = FALSE)
    mdata[,"sp_num"]<-as.factor(as.numeric(mdata[,"sp_num"]))
    mdata<-data.frame(mdata,stringsAsFactors = TRUE)
    
    if(pdf){
      # prepare plot size
      plen<-NULL
      if(is.null(pdflen)){
        # here need from data
        plen <- unique(sp_num)
        plen <- length(plen)
        plen <- plen * 0.7
        if(plen < 3){plen <- 3}
      } else {
        plen <- pdflen
      }
      pwid<-NULL
      if(is.null(pdfwid)){
        # here need from data
        pwid <- unique(sp)
        pwid <- length(pwid)
        pwid <- pwid * 0.7
        if(pwid < 3){pwid <- 3}
      } else {
        pwid <- pdfwid
      }
      # set self-defined colors 
      us<-unique(size)
      us<-us[order(us)]
      n<-length(us)
      tcol<-grDevices::cm.colors(n)
      cols<-rep(NA,length(size))
      for(ii in 1:length(us)){
        w<-which(size == us[ii])
        cols[w]<-tcol[ii]
      }
      cols
      
      # to make x-y dotplot
      g1<-
        ggplot2::ggplot(data = mdata, 
                        mapping = ggplot2::aes(x = sp,
                                               y = sp_num,
                                               size = size, 
                                               color = -size)) +
        ggplot2::geom_point() + 
        ggplot2::scale_color_distiller(name = "count", 
                              palette = 14) +
#                              breaks =levels(mdata$size)) +
#        ggplot2::geom_point(color = cols) + 
#        ggplot2::scale_color_brewer(pallete = col) +
        ggplot2::scale_size_area(max_size = 20, guide = FALSE) +
        ggplot2::theme(axis.text.x = element_text(angle = 90)) +
        ggplot2::scale_x_discrete(position="top") +
        ggplot2::scale_y_discrete(limits = rev(levels(mdata$sp_num))) +
        ggplot2::geom_text(ggplot2::aes(label = size), 
                           vjust = .5, colour = 'black', size = 3)
#      print(g1)
      pdfname1 <- paste("xyplot_of_species_distri_", 
                        pname, ".pdf",sep="")
      
      grDevices::pdf(pdfname1,pwid,plen)
      print(g1)
      grDevices::dev.off()
    }

    rr<-list(mdata)
    names(rr)<-"two-dimension"
    ress<-append(ress,rr)
  }

  
  #######################################################
  ## one-dimension plots
  #######################################################
  if(plot.type == "sp"){
    # make plot-1: x == each sp, y == count of gd
    allsp <- dat[ , "sp"]
    allsp <- unlist(strsplit(allsp, split = ",", 
                             fixed = TRUE))
    sp <- allsp[order(allsp)]
    data <- data.frame(sp)
    
    # make plot
    gglist <- vector("list", 2)
    g1 <- 
      ggplot2::ggplot(
        data, ggplot2::aes_string(
          x = stats::reorder(x = data[ , "sp"], 
                             X = rep(1, nrow(data)), 
                             FUN = function(x) sum(x))))
    
    gg1 <- g1 + ggplot2::geom_bar(ggplot2::aes_string(fill = "sp"),
                                  position = ggplot2::position_stack(reverse = TRUE),
                                  show.legend = FALSE) +
      ggplot2::coord_flip() + 
      ggplot2::guides(fill=ggplot2::guide_legend(title = "")) + 
      # no title for color legend
      ggplot2::xlab(label = "sp")
    

    if(pdf){
      plen<-NULL
      if(is.null(pdflen)){
        # here need from data
        plen <- unique(data[ ,"sp"])
        plen <- length(plen)
        plen <- plen * 0.2
        if(plen < 3){plen <- 3}
      } else {
        plen <- pdflen
      }
      pwid<-NULL
      if(is.null(pdfwid)){
        # here need from data
        pwid <- 5
      } else {
        pwid <- pdfwid
      }
      
      
      pdfname1 <- paste("Distribution_of_species_", 
                        pname, ".pdf",sep="")
      grDevices::pdf(file = pdfname1, width = pwid, 
                     height = plen)
      print(gg1)
      grDevices::dev.off()
    }
    
    rr<-list(data)
    names(rr)<-"sp"
    ress<-append(ress,rr)
  }
 
  
  
  if(plot.type == "sp_num"){
    # make plot-2: x == sp number, y == count of gd
    sp_num <- as.numeric(dat[ , "sp_num"])
    # make factor for order
    sp_num <- as.factor(as.numeric(sp_num))
    data <- data.frame(sp_num)
    
    # make plot
    g1 <- ggplot2::ggplot(data, ggplot2::aes_string(x = "sp_num"))
    
    gg1 <- g1 + ggplot2::geom_bar(
      ggplot2::aes_string(fill = "sp_num"), 
      position = ggplot2::position_stack(reverse = TRUE),
      show.legend = FALSE) +
      # make vertical plot
      ggplot2::coord_flip() + 
      # reverse order
      ggplot2::scale_x_discrete(limits = rev(levels(data$sp_num))) +
      ggplot2::guides(fill = ggplot2::guide_legend(title = "")) + 
      # no title for color legend
      ggplot2::xlab(label = "sp_num")
    
    gglist[[2]] <- gg1
    
    if(pdf){
      plen <- NULL
      if(is.null(pdflen)){
        # here grep from dat
        plen <- unique(dat[ ,"sp_num"])
        plen <- length(plen)
        plen <- plen * 0.2
        if(plen < 3){plen <- 3}
      } else {
        plen <- pdflen
      }
      pwid<-NULL
      if(is.null(pdfwid)){
        # here need from data
        pwid <- 5
      } else {
        pwid <- pdfwid
      }
      
      pdfname2 <- paste("Distribution_species_num_", 
                        pname, ".pdf", sep="")
      grDevices::pdf(file = pdfname2, width = pwid, 
                     height = plen)
      print(gg1)
      grDevices::dev.off()
    }
    
    rr<-list(data)
    names(rr)<-"sp_num"
    ress<-append(ress,rr)
  }
  
  return(ress)
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


