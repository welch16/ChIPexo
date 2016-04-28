

rm(list = ls())

library(data.table)
library(ggplot2)
library(scales)

data_dir <- "data/SCC_curves/"
files <- list.files(data_dir,include.dirs = TRUE,full.names = TRUE)

curves <- lapply(files,read.table,header = TRUE)
curves <- lapply(curves,data.table)

names(curves) <- gsub("_SCC.txt","",basename(files))

rangeT <- mapply(function(curve,name){
  rr <- curve[,range(cross.corr)]
  out <- data.table(name , min = rr[1],max = rr[2])
  out
  },curves,names(curves),SIMPLIFY = FALSE)

rangeT <- do.call(rbind,rangeT)
rangeT[,diff := max - min]
rangeT[,nsc := max / min]


plot_curve <- function(x,name){
  ggplot(x,aes(shift,cross.corr))+geom_line()+
    ggtitle(name)+scale_x_continuous(breaks = seq(10,300,by = 20))+
    xlab("Shift")+ylab("Strand Cross - Correlation (SCC)")
}

plots <- mapply(plot_curve,curves,names(curves),SIMPLIFY = FALSE)

figs_dir <- "figs/all_SCC"


save_plot <- function(pp,file){
  pdf(file,width = 6 , height  = 4)
  print(pp)
  dev.off()
}

outfiles <- file.path(figs_dir,gsub(".txt",".pdf",basename(files)))

mapply(save_plot,plots,outfiles)
