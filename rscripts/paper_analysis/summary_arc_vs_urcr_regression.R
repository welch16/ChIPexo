
rm(list = ls())

library(broom)
library(data.table)
library(ggplot2)
library(scales)
library(RColorBrewer)

in_dir <- "data/ARC_quantify"
files <- list.files(in_dir)

load_file <- function(x){load(x);coeff}

coeffs <- lapply(file.path(in_dir,files),load_file)

coeffs <- mapply(function(x,y){
  x$id = y
  x},coeffs,gsub(".RData","",files),SIMPLIFY = FALSE)


coeffs <- data.table(do.call(rbind,coeffs))


coeffs[ , Rep := "Rep-1"]

coeffs[grep("Rep2",id), Rep := "Rep-2"]
coeffs[grep("Rep3",id), Rep := "Rep-3"]

coeffs[grep("1314",id), Rep := "Rep-2"]
coeffs[grep("1320",id), Rep := "Rep-2"]
coeffs[grep("933",id), Rep := "Rep-2"]
coeffs[grep("937",id), Rep := "Rep-2"]

coeffs[, expt := "a"]
coeffs[ Rep == "Rep-1", expt := gsub("_Rep1","",id)]
coeffs[ Rep == "Rep-2", expt := gsub("_Rep2","",id)]
coeffs[ Rep == "Rep-3", expt := gsub("_Rep3","",id)]

coeffs[grepl("landick_sig70_9",id), expt := "landick_aero"]
coeffs[grepl("landick_sig70_13",id), expt := "landick_sig70_rif"]

coeffs[,organism := "Human"]
coeffs[grepl("mouse",id),organism := "Mouse"]
coeffs[grepl("landick",id),organism := "E.Coli"]
coeffs[grepl("S2",id),organism := "D.Melanogaster"]
coeffs[grepl("embryo",id),organism := "D.Melanogaster"]

coeffs[ , seq := "ChIP-exo"]
coeffs[ grepl("nexus",id), seq := "ChIP-nexus"]
coeffs <- coeffs[order(seq,organism)]


z <- qnorm(.95)

figs_dir <- "figs/test_quantify"

pdf(file = file.path(figs_dir,"arc_urcr_param_regression.pdf"),width = 12,height = 6)
ggplot(coeffs[term == "b"],aes(id,estimate,fill = organism))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(angle = 0),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0))+
  facet_grid( ~ seq + organism  ,scales = "free",space = "free")+
  ylab("Intercept")+xlab("Experiment")+
  scale_fill_brewer(palette = "Set1")+ggtitle("B")
ggplot(coeffs[term == "k"],aes(id,estimate,fill = organism))+
  geom_bar(stat = "identity")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(angle = 0),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0))+
  facet_grid( ~ seq + organism  ,scales = "free",space = "free")+
  ylab("Slope")+xlab("Experiment")+
  scale_fill_brewer(palette = "Set1")+ggtitle("C")
dev.off()

  ## geom_text(data=coeffs[term == "b"],aes(id,estimate,label=round(estimate,2)),vjust=1,hjust = .5,angle = 90) +

  ## geom_text(data=coeffs[term == "k"],aes(id,estimate,label=round(estimate,2)),vjust=1,hjust = .5,angle = 90) +
