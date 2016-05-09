
rm(list = ls())

library(scales)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(devtools)
library(viridis)
library(gridExtra)

data_dir <- "data/npos_width_reg_sig70"
files <- list.files(data_dir,pattern = "RData",full.names = TRUE,
                    include.dirs = TRUE)
nms <- gsub(".RData","",basename(files))

load_coeff <- function(x){
  load(x)
  return(coeff)}

coeffs <- lapply(files,load_coeff)
names(coeffs) <- nms

coeffs <- mapply(function(x,y)x[,name := y],
         coeffs,names(coeffs),SIMPLIFY = FALSE)
coeffs <- do.call(rbind,coeffs)

meta <- data.table(name = nms)


get_sample <- function(name){
  name <- strsplit(name,"_")
  name <- lapply(name,function(x){
    qq <- grepl("subsam",x)
    if(any(qq)){
      out <- gsub("subsample","",x[qq])
      out <- as.numeric(out)
    }else{
      out <- "complete"
    }
    return(out)
  })
  return(do.call(c,name))
}

get_edsn <- function(name){
  name <- strsplit(name,"_")
  name <- lapply(name,function(x){
    qq <- grepl("edsn",x)
    if(any(qq)){
      out <- x[qq]
    }else{
      out <- ""
    }
    return(out)})
  return(do.call(c,name))
}
     

meta <- meta %>% mutate(
   samp = get_sample(name),
   edsn = get_edsn(name))


## K562 - TBP

figs_dir <- "figs/for_paper"


coeffs <- merge(coeffs,meta , by = "name")

r <- brewer.pal(8,name = "Dark2")

DT <- coeffs


pdf(file.path(figs_dir,"sig70_qc_boxplot_eval.pdf"),width = 16,height = 4)
p1 <- ggplot(DT[term == "npos"],aes(as.factor(samp),estimate,
                  colour = edsn))+
  geom_boxplot(outlier.size = NA)+
  facet_grid(. ~ edsn ,scales = "free",space = "free")+  
  scale_y_log10()+ylab("npos")+
  theme_bw()+theme(legend.position = "none",
                   axis.text.x = element_text(angle = 90),
                   axis.title.x = element_blank(),
                   plot.title = element_text(hjust = 0),
                   strip.background = element_blank())+
  geom_abline(slope = 0 ,intercept = 1,linetype = 2 ,colour = "darkgrey")+
  ylab("Adjusted Average Read Coverage")
p2 <- ggplot(DT[term == "width"],aes(as.factor(samp),-estimate,colour = edsn))+
  geom_boxplot(outlier.size = NA)+
  facet_grid(.~  edsn ,scales = "free",space = "free")+
  theme_bw()+theme(legend.position = "none",
                   axis.text.x = element_text(angle = 90),
                   axis.title.x = element_blank(),
                   plot.title = element_text(hjust = 0),
                   strip.background = element_blank())+
  ylim(-.1,.25)+ylab("Average Read Coverage Bias")+
  geom_abline(slope = 0 ,intercept = 0,linetype = 2 ,colour = "darkgrey")
p1
p2
dev.off()


pdf(file.path(figs_dir,"sig70_qc_line_eval.pdf"),height = 5,width = 5)
DT1 <- DT[ ,mean(estimate), by = .(term,samp,edsn)]
DT2 <- DT[ ,median(estimate), by = .(term,samp,edsn)]
p1 <- ggplot(DT1[term == "npos"],aes(samp,V1,colour = edsn))+
  geom_line()+
  theme_bw()+theme(axis.text.x = element_text(angle = 90),
                   axis.title.x = element_blank(),
                   plot.title = element_text(hjust = 0),
                   strip.background = element_blank())+
  geom_abline(slope = 0 ,intercept = 1,linetype = 2 ,colour = "darkgrey")+
  ylab("Adjusted Average Read Coverage")+
  scale_color_brewer(palette = "Dark2")
p2 <- ggplot(DT1[term == "width"],aes(samp,-V1,colour = edsn))+
  geom_line()+
  theme_bw()+theme(axis.text.x = element_text(angle = 90),
                   axis.title.x = element_blank(),
                   plot.title = element_text(hjust = 0),
                   strip.background = element_blank())+
  ylim(0,.1)+ylab("Average Read Coverage Bias")+
  geom_abline(slope = 0 ,intercept = 0,linetype = 2 ,colour = "darkgrey")+
  scale_color_brewer(palette = "Dark2")
p1
p2
dev.off()

