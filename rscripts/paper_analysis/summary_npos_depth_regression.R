
rm(list = ls())

library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)

data_dir <- "data/QC_quantify"
files <- list.files(data_dir)

load_file <- function(x){
  load(x)
  return(coeff)
}

coeffs <- lapply(file.path(data_dir,files),load_file)
coeffs <- mapply(function(x,y)x[,id := y],coeffs,gsub(".RData","",files),SIMPLIFY = FALSE)
DT <- do.call(rbind,coeffs)


N <- 1e4
K <- N

mean_est <- sapply(coeffs,function(x)x[,mean(estimate[1:K])])
var_est <- sapply(coeffs,function(x)x[,var(estimate[1:K])])

z <- qnorm(.05,lower.tail = FALSE)

dt <- data.table(id = sub(".RData","",files),
                 m = mean_est,var = var_est)
dt[,lb := m - z * sqrt(var) / N]
dt[,ub := m + z * sqrt(var) / N]
dt[,organism := "Human"]
dt[grepl("mouse",id),organism := "Mouse"]
dt[grepl("landick",id),organism := "E.Coli"]
dt[grepl("S2",id),organism := "D.melanogaster"]
dt[grepl("embryo",id),organism := "D.melanogaster"]

DT[,organism := "Human"]
DT[grepl("mouse",id),organism := "Mouse"]
DT[grepl("landick",id),organism := "E.Coli"]
DT[grepl("S2",id),organism := "D.Melanogaster"]
DT[grepl("embryo",id),organism := "D.Melanogaster"]





figs_dir <- "figs/test_quantify"
pdf(file = file.path(figs_dir,"slope_npos_depth.pdf"),width = 12,height = 6)
ggplot(dt,aes(id,m))+geom_bar(stat = "identity")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")
ggplot(dt[!grepl("venter",id)],aes(id,m))+geom_bar(stat = "identity")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")
ggplot(dt[!grepl("venter",id)],aes(id,m))+geom_bar(stat = "identity")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")
ggplot(dt[grepl("carroll",id)],aes(id,m))+geom_bar(stat = "identity")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")
ggplot(dt[organism == "human"],aes(id,m))+geom_bar(stat = "identity")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")
ggplot(dt[grepl("landick",id)],aes(id,m))+geom_bar(stat = "identity")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")
ggplot(dt[grepl("nexus",id)],aes(id,m))+geom_bar(stat = "identity")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")
dev.off()

  ## 

DT[ , seq := "ChIP-exo"]
DT[ grepl("nexus",id), seq := "ChIP-nexus"]


DT <- DT[order(seq,organism)]


figs_dir <- "figs/test_quantify"
pdf(file = file.path(figs_dir,"slope_npos_depth_boxplot.pdf"),width = 12,height = 6)
ggplot(DT,aes(id,estimate,fill = organism , colour = organism))+
  geom_boxplot(outlier.size = 0)+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90),
        strip.text.x = element_text(angle = 0),
        strip.background = element_blank(),
        plot.title = element_text(hjust= 0))+
  facet_grid( ~ seq + organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Ratio between depth and nr. unique positions change")+
  scale_y_log10()+scale_fill_brewer(palette = "Pastel1")+scale_color_brewer(palette = "Set1")+
  geom_abline(slope = 0, intercept = 1,linetype = 2,colour = "darkgrey")+ggtitle("A")
ggplot(DT,aes(id,estimate))+geom_boxplot(outlier.size = 0)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")+
  scale_y_log10()
ggplot(DT[!grepl("venter",id)],aes(id,estimate))+geom_boxplot(outlier.size = 0)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")+
  scale_y_log10()
ggplot(DT[!grepl("venter",id)],aes(id,estimate))+geom_boxplot(outlier.size = 0)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")+
  scale_y_log10()
ggplot(DT[grepl("carroll",id)],aes(id,estimate))+geom_boxplot(outlier.size = 0)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")+
  scale_y_log10()
ggplot(DT[organism == "human"],aes(id,estimate))+geom_boxplot(outlier.size = 0)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")+
  scale_y_log10()
ggplot(DT[grepl("landick",id)],aes(id,estimate))+geom_boxplot(outlier.size = 0)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")+
  scale_y_log10()
ggplot(DT[grepl("nexus",id)],aes(id,estimate))+geom_boxplot(outlier.size = 0)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90))+
  facet_grid( ~ organism  ,scales = "free",space = "free")+xlab("")+
  xlab("Experiment")+ylab("Depth to Nr. of Unique positions rate")+
  scale_y_log10()
dev.off()


