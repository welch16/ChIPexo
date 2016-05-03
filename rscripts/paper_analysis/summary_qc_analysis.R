
rm(list = ls())

library(data.table)

data_dir <- "data/ChIPseq_QC"


genome_length <- data.table(genome = c("hg19","mm9","dm3","ecoli.K12"),
                            length = c(3095693981,2725765481,175e6,4639221))

g <- genome_length[,(length)]

genome_length[,ff := g[1]/ length]

## > g[1] / g[2] human to mouse
## [1] 1.135715 
## > g[1] / g[3] human to dm3
## [1] 17.68968
## > g[1] / g[4] human to ecoli
## [1] 667.2875


files <- list.files(data_dir,recursive = TRUE,full.names = TRUE,
                    include.dirs = TRUE)

qc_ind <- lapply(files,read.table,header = TRUE)
qc_ind <- lapply(qc_ind,data.table)

qc_ind <- mapply(function(x,y)x[,name := y],qc_ind,
                 gsub("_QC.txt","",basename(files)),
                 SIMPLIFY = FALSE)
qc_ind <- do.call(rbind,qc_ind)
qc_ind[,se := "exo"]
qc_ind[grep("nexus",name), se := "nexus"]
qc_ind[grep("SE",name), se := "chipseq(SE)"]
qc_ind[grep("PE",name), se := "chipseq(PE)"]

qc_ind[,genome := "hg19"]
qc_ind[grep("landick",name),genome := "ecoli.K12"]
qc_ind[grep("mouse",name),genome := "mm9"]
qc_ind[grep("embryo",name),genome := "dm3"]
qc_ind[grep("S2",name),genome := "dm3"]

qc_ind <- merge(qc_ind,genome_length,by = "genome")

qc_ind[,norm_depth := depth * ff]
qc_ind[,length := NULL]
qc_ind[,ff := NULL]


## library(ggplot2)
## library(ggrepel)

## ggplot(qc_ind,aes(norm_depth,NSC,colour = genome))+
##   geom_point()+
##   geom_label_repel(aes(label = name))+
##   scale_color_brewer(palette = "Set1")+scale_x_log10()+
##   theme(legend.position = "none")+ylim(0,75)
## dev.off()
