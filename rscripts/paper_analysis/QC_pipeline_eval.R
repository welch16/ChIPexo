
rm(list = ls())

library(scales)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(devtools)
library(viridis)
library(gridExtra)

data_dir <- "data/npos_width_reg"
files <- list.files(data_dir,pattern = "RData",full.names = TRUE,
                    include.dirs = TRUE)
nms <- gsub("_factors.RData","",basename(files))
nms <- gsub("mouse","FoxA1_liver",nms)
nms <- gsub("human","ER_MCF-7",nms)
nms <- gsub("CTCF","CTCF_HeLa_Rep1",nms)

load_coeff <- function(x){
  load(x)
  return(coeff)}

coeffs <- lapply(files,load_coeff)
names(coeffs) <- nms

coeffs <- mapply(function(x,y)x[,name := y],
         coeffs,names(coeffs),SIMPLIFY = FALSE)
coeffs <- do.call(rbind,coeffs)

meta <- data.table(name = nms)

get_repl <- function(name){
  repl <- strsplit(name,"_")
  repl <- lapply(repl,function(x){
    if(any(grepl("Rep",x))){
      out <- x[grepl("Rep",x)]
    }else{
      out <- "Rep1"
    }
    return(out) })
  repl <- do.call(c,repl)
  id1 <- grepl("landick",name) & grepl("sig70_9",name)
  repl[id1] <- paste0("Rep",rep(1:2,2))
  id2 <- grepl("landick",name) & grepl("sig70_1",name)
  repl[id2] <- paste0("Rep",rep(1:2,each = 2) )
  return(repl)
}

get_genome <- function(name)ifelse(grepl("landick",name),"ecoli.K12",
   ifelse(grepl("liver",name),"mm9",
   ifelse(grepl("S2",name),"dm3",
   ifelse(grepl("embryo",name),"dm3","hg19"))))

get_sample <- function(name){
  name <- strsplit(name,"_")
  name <- lapply(name,function(x){
    qq <- grepl("subsam",x)
    if(any(qq)){
      out <-paste0( gsub("subsample","",x[qq]),"M")
    }else{
      out <- "complete"
    }
    return(out)
  })
  return(do.call(c,name))
}

get_author <- function(name){
  ifelse(grepl("carroll",name),"Carroll",
  ifelse(grepl("nexus",name),"Zeitlinger",
  ifelse(grepl("venters",name),"Pugh",
  ifelse(grepl("land",name),"Landick",
  ifelse(grepl("pugh",name),"Pugh","Meijsing")))))}

get_TF <- function(name){
  out <- ifelse(grepl("carroll_hu",name),"ER",
  ifelse(grepl("carroll_mo",name),"FoxA1",
  ifelse(grepl("land",name),"Sig70",
  ifelse(grepl("vent",name),"TBP",
  ifelse(grepl("CTCF",name),"CTCF","")))))

  id1 <- grepl("meij",name)
  out[id1] <- sapply(strsplit(name[id1],"_"),function(x)x[2])

  id2 <- grepl("nexus",name)
  out[id2] <- sapply(strsplit(name[id2],"_"),function(x)x[3])

  
  return(out)
}

get_cell <- function(name){
  out <- ifelse(grepl("carroll_hu",name),"MCF-7",
  ifelse(grepl("carroll_mo",name),"liver",
  ifelse(grepl("land",name),"-",
  ifelse(grepl("vent",name),"K562",
  ifelse(grepl("CTCF",name),"HeLA","")))))

  id1 <- grepl("meij",name)
  out[id1] <- sapply(strsplit(name[id1],"_"),function(x)x[3])

  id2 <- grepl("nexus",name)
  out[id2] <- sapply(strsplit(name[id2],"_"),function(x)x[2])

  return(out)
}



meta <- meta %>% mutate(
   protocol = ifelse(grepl("nexus",name),"ChIP-nexus","ChIP-exo"),
   genome = get_genome(name),
   repl = get_repl(name),
   samp = get_sample(name),
   lab = get_author(name),
   TF = get_TF(name),
   cell = get_cell(name))

meta[, name2 := sapply(strsplit(meta$name,"_"),function(x)do.call(paste,as.list(x[-1])))]


## K562 - TBP

DT <- coeffs[grepl("K562",name)][!grepl("mei",name)]


DT <- merge(DT,meta,by = "name") %>%
    filter(name != "chipnexus_K562_TBP_Rep1")


  ## filter(!grepl("venters_TBP_K562_Rep1",name))



figs_dir <- "figs/for_paper"


coeffs <- merge(coeffs,meta , by = "name")

by_term_name <-
  group_by(coeffs,term,name2)
qc <- summarize(by_term_name,
                mean = mean(estimate),
                median = median(estimate),
                sd = sd(estimate),
                gamma = 1 / mean(estimate))


z <- qnorm(.95)



DT2 <- coeffs[samp == "complete"][lab != "Landick"]

r <- brewer.pal(8,name = "Set1")

pdf(file = file.path(figs_dir,"QC_pipeline_eval_boxplot.pdf"),width = 12,height = 6)
ggplot(DT2[term == "npos"],aes(name2, estimate,colour = interaction(genome,lab) ))+
  geom_boxplot(outlier.size = NA)+
  facet_grid(.  ~genome + protocol +lab,scales = "free",space = "free")+theme_bw()+
  scale_y_log10()+geom_abline(slope = 0,intercept = 1,linetype = 2,colour = "darkgrey")+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0))+
  scale_color_manual(values = r[-6])+ggtitle("A")+
  ylab("Adjusted Average Read Coverage")
ggplot(DT2[term == "width"],aes(name2, -estimate,
             colour = interaction(genome,lab)))+
  geom_boxplot(outlier.size = NA)+
  facet_grid(. ~ genome + protocol +lab,scales = "free",space = "free")+theme_bw()+
  geom_abline(slope = 0,intercept = 0,linetype = 2,colour = "darkgrey")+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust= 0))+ggtitle("B")+
  scale_color_manual(values = r[-6])+
  ylab("Average Read Coverage Bias")+ylim(-1,5)
dev.off()



pdf(file = file.path(figs_dir,"QC_pipeline_eval_boxplot2.pdf"),width = 12,height = 5)
ggplot(DT2[term == "npos"],aes(name2,1/ estimate,colour = interaction(genome,lab) ))+
  geom_boxplot(outlier.size = NA)+
  facet_grid(.  ~genome + protocol +lab,scales = "free",space = "free")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = "none",
        strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0))+
  scale_color_manual(values = r[-6])+ggtitle("A")+ylab("1/beta_1")+ylim(0,1)
dev.off()





pdf(file = file.path(figs_dir,"K562_TBP_sample.pdf"),width = 6,height = 4)
p1 <- ggplot(DT[term == "npos"],aes(samp,estimate,
                  colour = paste0(protocol,repl)))+
  geom_boxplot(outlier.size = NA)+
  facet_grid( . ~  protocol + repl  )+
  scale_y_log10()+ylab("npos")+
  theme_bw()+theme(legend.position = "none",
                   axis.text.x = element_text(angle = 90),
                   axis.title.x = element_blank(),
                   plot.title = element_text(hjust = 0),
                   strip.background = element_blank())+
  scale_color_manual(values = r[c(4,4,4,7)])+ggtitle("C")+  
  geom_abline(slope = 0 ,intercept = 1,linetype = 2 ,colour = "darkgrey")+
  ylab("Adjusted Average Read Coverage")
p2 <- ggplot(DT[term == "width"],aes(samp,-estimate,
                  colour = paste0(protocol,repl)))+
  geom_boxplot(outlier.size = NA)+
  facet_grid( . ~  protocol + repl  )+ylab("width")+
  theme_bw()+theme(legend.position = "none",
                   axis.text.x = element_text(angle = 90),
                   axis.title.x = element_blank(),
                   plot.title = element_text(hjust = 0),
                   strip.background = element_blank())+
  ylim(-1,5)+ylab("Average Read Coverage Bias")+
  scale_color_manual(values = r[c(4,4,4,7)])+ggtitle("D")+
  geom_abline(slope = 0 ,intercept = 0,linetype = 2 ,colour = "darkgrey")
p1
p2
dev.off()

DT <- DT %>% filter(samp != "complete")


DT <- DT %>% select(term,estimate,protocol,repl,samp,lab)
coeff <- group_by(DT,term,protocol,repl,samp)
out <- summarize(coeff,
   med = ifelse(term == "npos",1,-1) * median(estimate),
   mean = ifelse(term == "npos",1,-1) * mean(estimate),
   sd = sd(estimate),
   trim = ifelse(term == "npos",1,-1) * mean(estimate,trim = .1))

   ## , 


out <- out %>% mutate(samp2 = as.numeric(gsub("M","",samp)))

library(gridExtra)


p1 <- ggplot(out[term == "npos"],aes(samp2,
             mean,
             shape = repl,
             linetype = repl,
             colour = protocol))+
  geom_point(size = 3)+geom_line(size = .5)+scale_y_log10()+
  scale_linetype_discrete(guide = "none")+
  geom_abline(slope = 0,intercept = 1,linetype = 2)+
  scale_color_manual(values = r[c(4,7)],guide = "none")+
  scale_shape_discrete(solid = FALSE,name = "Replicate")+
  theme_bw()+theme(legend.position = "top")+
  scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
  xlab("Number of reads")+ylab("param1")                              
p2 <- ggplot(out[term == "width"],aes(samp2,
             mean,
             shape = repl,
             linetype = repl,
             colour = protocol))+
  geom_point(size = 3)+geom_line(size = .5)+
  scale_linetype_discrete(guide = "none")+
  geom_abline(slope = 0,intercept = 0,linetype = 2)+
  scale_color_manual(values = r[c(4,7)],guide = "none")+
  scale_shape_discrete(solid = FALSE,name = "Replicate")+
  theme_bw()+theme(legend.position = "top")+
  scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
  xlab("Number of reads")+ylab("param2")
m1 <- ggplot(out[term == "npos"],aes(samp2,
             med,
             shape = repl,
             linetype = repl,
             colour = protocol))+
  geom_point(size = 3)+geom_line(size = .5)+scale_y_log10()+
  scale_linetype_discrete(guide = "none")+
  geom_abline(slope = 0,intercept = 1,linetype = 2)+
  scale_color_manual(values = r[c(4,7)],guide = "none")+
  scale_shape_discrete(solid = FALSE,name = "Replicate")+
  theme_bw()+theme(legend.position = "top")+
  scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
  xlab("Number of reads")+ylab("param1")
m2 <- ggplot(out[term == "width"],aes(samp2,
             med,
             shape = repl,
             linetype = repl,
             colour = protocol))+
  geom_point(size = 3)+geom_line(size = .5)+
  scale_linetype_discrete(guide = "none")+
  geom_abline(slope = 0,intercept = 0,linetype = 2)+
  scale_color_manual(values = r[c(4,7)],guide = "none")+
  scale_shape_discrete(solid = FALSE,name = "Replicate")+
  theme_bw()+theme(legend.position = "top")+
  scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
  xlab("Number of reads")+ylab("param2")
t1 <- ggplot(out[term == "npos"],aes(samp2,
             trim,
             shape = repl,
             linetype = repl,
             colour = protocol))+
  geom_point(size = 3)+geom_line(size = .5)+scale_y_log10()+
  scale_linetype_discrete(guide = "none")+
  geom_abline(slope = 0,intercept = 1,linetype = 2)+
  scale_color_manual(values = r[c(4,7)],guide = "none")+
  scale_shape_discrete(solid = FALSE,name = "Replicate")+
  theme_bw()+theme(legend.position = "top")+
  scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
  xlab("Number of reads")+ylab("param1")
t2 <- ggplot(out[term == "width"],aes(samp2,
             trim,
             shape = repl,
             linetype = repl,
             colour = protocol))+
  geom_point(size = 3)+geom_line(size = .5)+
  scale_linetype_discrete(guide = "none")+
  geom_abline(slope = 0,intercept = 0,linetype = 2)+
  scale_color_manual(values = r[c(4,7)],guide = "none")+
  scale_shape_discrete(solid = FALSE,name = "Replicate")+
  theme_bw()+theme(legend.position = "top")+
  scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
  xlab("Number of reads")+ylab("param2")
q1 <- ggplot(out[term == "npos"],aes(samp2,
             sd,
             shape = repl,
             linetype = repl,
             colour = protocol))+
  geom_point(size = 3)+geom_line(size = .5)+scale_y_log10()+
  scale_linetype_discrete(guide = "none")+
  geom_abline(slope = 0,intercept = 1,linetype = 2)+
  scale_color_manual(values = r[c(4,7)],guide = "none")+
  scale_shape_discrete(solid = FALSE,name = "Replicate")+
  theme_bw()+theme(legend.position = "top")+
  scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
  xlab("Number of reads")+ylab("param1")
q2 <- ggplot(out[term == "width"],aes(samp2,
             sd,
             shape = repl,
             linetype = repl,
             colour = protocol))+
  geom_point(size = 3)+geom_line(size = .5)+
  scale_linetype_discrete(guide = "none")+
  geom_abline(slope = 0,intercept = 0,linetype = 2)+
  scale_color_manual(values = r[c(4,7)],guide = "none")+
  scale_shape_discrete(solid = FALSE,name = "Replicate")+
  theme_bw()+theme(legend.position = "top")+
  scale_x_continuous(labels = paste0(seq(20,50,by = 10),"M"))+
  xlab("Number of reads")+ylab("param2")


pdf(file = file.path(figs_dir,"TBP_param_depth_trend.pdf"),width = 8,height=4)
grid.arrange(p1,p2,nrow = 1)
grid.arrange(m1,m2,nrow = 1)
grid.arrange(t1,t2,nrow = 1)
grid.arrange(q1,q2,nrow = 1)
dev.off()

pdf(file = file.path(figs_dir,"TBP_param_depth_trend2.pdf"),width = 4,height=4)
p1
p2
m1
m2
t1
t2
q1
q2
dev.off()
