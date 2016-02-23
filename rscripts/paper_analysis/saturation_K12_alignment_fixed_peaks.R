
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ChIPUtils)
library(RColorBrewer)
library(ggplot2)
library(scales)

indir <- "/p/keles/ChIPexo/volume6/K12/saturation"
figsdir <- "figs/saturation/K12_alignment"

all_files <- list.files(indir,recursive = TRUE)
all_files <- all_files[grep("bam",all_files,invert = TRUE)]

annot_dir <- "/p/keles/ChIPexo/volume6/K12/annotations"
afiles <- list.files(annot_dir)
afiles <- afiles[grep("csv",afiles)]

annots <- read.csv(file.path(annot_dir,afiles),sep = "\t")
annots <- data.table(annots)
annots <- annots[grep("Sigma70",sigma.factor)]
annots <- annots[grep("TIM",description)]
annots <- annots[grep("RPP",description,invert = TRUE)]
annots <- annots[grep("AIPP",description,invert = TRUE)]
annots <- annots[grep("WHO",description,invert = TRUE)]
annots <- annots[grep("NTAS",description,invert = TRUE)]
annots <- annots[grep("ICA",description,invert = TRUE)]
annotDT <- annots

annots <- annots[,IRanges(start = coord,width = 1)]
annots <- sort(annots)


### test example
edsn <- data.table(exo = c(1311,1314,1317,1320),seq = c(1396,1398,1400,1402))
seeds <- list.files(file.path(indir,"ChIPexo"))
seeds <- gsub("seed","",seeds)

samps <- seq(100,900,by = 100)
ext <- 15
max_dist <- 100
mm <- 300

center <- mid(annots)
annots <- IRanges(center - ext ,center + ext)

uni_peaks <- function(files){
  out <- lapply(files,read.table)
  out <- lapply(out,data.table)
  out <- lapply(out,function(x)x[,IRanges(start = V2,end = V3)])
  out <- sort(do.call(intersect,out))
  return(out)
}
  


filter_annots <- function(annots,edsn,edsn_idx,exo_files,pet_files,set_files)
{
  exo_peaks <- exo_files[grepl("peak",exo_files) & grepl("900K",exo_files) &
                         grepl(edsn[edsn_idx,(exo)],exo_files)]
  pet_peaks <- pet_files[grepl("peak",pet_files) & grepl("900K",pet_files) &
                         grepl(edsn[edsn_idx,(seq)],pet_files)]
  set_peaks <- set_files[grepl("peak",set_files) & grepl("900K",set_files) &
                         grepl(edsn[edsn_idx,(seq)],set_files)]

  exo_peaks <- uni_peaks(file.path(indir,exo_peaks))
  pet_peaks <- uni_peaks(file.path(indir,pet_peaks))  
  set_peaks <- uni_peaks(file.path(indir,set_peaks))

  peaks <- exo_peaks[ countOverlaps(exo_peaks,pet_peaks) > 0 & countOverlaps(exo_peaks,set_peaks) > 0]

  annots <- subsetByOverlaps(annots,peaks)

  return(annots)
}

simple_saturation_analysis <- function(samp,edsn_idx,seed,edsn,all_files,annots,topM = 500)
{
  exo_files <- all_files[grep("ChIPexo",all_files)]
  pet_files <- all_files[grep("ChIPseq_PET",all_files)]
  set_files <- all_files[grep("ChIPseq_SET",all_files)]

   ## filter annots
  annots <- filter_annots(annots,edsn,edsn_idx,exo_files,pet_files,set_files)

  samptxt <- paste0("samp",samp)

  peakfiles <- list()
  peakfiles[["exo"]] <- exo_files[grepl("900K",exo_files) & grepl("peak",exo_files) &
                                  grepl(edsn[edsn_idx,(exo)],exo_files) & grepl(seed,exo_files)]
  peakfiles[["pet"]] <- pet_files[grepl("900K",pet_files) & grepl("peak",pet_files) &
                                  grepl(edsn[edsn_idx,(seq)],pet_files) & grepl(seed,pet_files)]
  peakfiles[["set"]] <- set_files[grepl("900K",set_files) & grepl("peak",set_files) &
                                  grepl(edsn[edsn_idx,(seq)],set_files) & grepl(seed,set_files)]

  exo_files <- exo_files[grepl(samptxt,exo_files) & grepl(edsn[edsn_idx,(exo)],exo_files) &
                         grepl(seed,exo_files) & grepl("fixed",exo_files)]
  pet_files <- pet_files[grepl(samptxt,pet_files) & grepl(edsn[edsn_idx,(seq)],pet_files) &
                         grepl(seed,pet_files) & grepl("fixed",pet_files)]
  set_files <- set_files[grepl(samptxt,set_files) & grepl(edsn[edsn_idx,(seq)],set_files) &
                         grepl(seed,set_files) & grepl("fixed",set_files)]

  peaks <- lapply(file.path(indir,peakfiles),read.table)
  peaks <- lapply(peaks,data.table)
  names(peaks) <- names(peakfiles)
  
  sites <- list()
  sites[["exo"]] <- data.table(read.table(file.path(indir,exo_files),header = TRUE))
  sites[["pet"]] <- data.table(read.table(file.path(indir,pet_files),header = TRUE))
  sites[["set"]] <- data.table(read.table(file.path(indir,set_files),header = TRUE))

  peaks <- lapply(peaks,function(x)x[order(-V8)][1:ifelse(topM > nrow(x),nrow(x),topM) ])
  peak_ranges <- lapply(peaks,function(x)x[,IRanges(start = V2,end = V3)])
  site_ranges <- lapply(sites,function(x)x[,IRanges(start = start,end = end)])
  overlaps <- mapply(findOverlaps,peak_ranges,site_ranges,SIMPLIFY = FALSE)
  site_ranges <- mapply(subsetByOverlaps,site_ranges,peak_ranges,SIMPLIFY = FALSE)
  nPeak <- sapply(overlaps,function(x)length(unique(queryHits(x))))
  overlaps <- mapply(findOverlaps,peak_ranges,site_ranges,SIMPLIFY = FALSE)
  nSite <- mapply(function(site,ov)length(site[subjectHits(ov)]),site_ranges,overlaps)

  nIden <- sapply(site_ranges,function(x,ann)length(findOverlaps(x,ann)),annots)
  annots <- subsetByOverlaps(annots,peak_ranges[["exo"]])
  reso <- sapply(site_ranges,function(x,ann){
    reso <- sapply(ann,function(y){
      d <- min(abs(y - mid(x)))
      if(d > max_dist){
        d <- NA
      }
      return(d)})
    median(reso,na.rm = TRUE)},mid(annots))
  out <- data.table(seq = names(peaks),nPeak,nSite,nIden,reso, N = samp * 1e3,condition = edsn_idx,seed)
  return(out)
}


param <- data.table(expand.grid(k = 1:4,samp = samps,seed = seeds))
saturation <- mclapply(1:nrow(param),
  function(i){
    sat <- simple_saturation_analysis(param[i,(samp)],param[i,(k)],param[i,(seed)],
         edsn,all_files,annots,topM = mm)
    if(!is.na(sat) && nrow(sat) > 1 ){
      return(sat)      
    }else{
      return(NULL)
    }},mc.cores = detectCores())
    

saturation <- do.call(rbind,saturation)

sat <- saturation[,{
  nPeak = round(mean(nPeak),0)
  nSite = round(mean(nSite),0)
  nIden = round(mean(nIden),0)
  resol = round(mean(reso),0)
  list(nPeaks = nPeak,nSites = nSite,nIden = nIden,reso = resol)},by = .(seq,N,condition)]

sat[, seq := plyr::mapvalues(seq, from = c("exo","pet","set"),
        to = c("ChIP-exo","PE ChIP-Seq","SE ChIP-Seq"))]

sat[,rif := ifelse(condition %in% c(1,3),"0min","20min")]
sat[,repl := ifelse(condition %in% 1:2,"Rep-1","Rep-2")]

library(MASS)
library(splines)

pdf(file = file.path(figsdir,"saturation_plots_panel_peaksFixed.pdf"))
ggplot(sat,aes(N,nPeaks,colour = seq))+geom_line()+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  facet_grid( repl ~ rif)+theme_bw()+
  theme(legend.position ="top")+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of candidate regions")+xlab("Nr. of reads") 
ggplot(sat,aes(N,nSites,colour = seq))+stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3))+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  facet_grid( repl ~ rif)+theme_bw()+
  theme(legend.position ="top")+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of predicted events")+xlab("Nr. of reads")
ggplot(sat,aes(N,nIden,colour = seq))+stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3))+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  facet_grid( repl ~ rif)+theme_bw()+
  theme(legend.position ="top")+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of identified targets")+xlab("Nr. of reads")
ggplot(sat,aes(N,reso,colour = seq))+stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3))+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  facet_grid( repl ~ rif)+theme_bw()+
  theme(legend.position ="top")+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Resolution")+xlab("Nr. of reads")
dev.off()


pdf(file = file.path(figsdir,"Sig70_rep1_rif0_saturation_analysis.pdf"))
ggplot(sat[ rif == "0min" & repl == "Rep-1"],aes(N,nPeaks,colour = seq))+
  geom_line(size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of candidate regions")+xlab("Nr. of reads") 
ggplot(sat[ rif == "0min" & repl == "Rep-1"],aes(N,nSites,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylim(380,560)+
  ylab("Number of predicted events")+xlab("Nr. of reads")
ggplot(sat[ rif == "0min" & repl == "Rep-1"],aes(N,nIden,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(55,120)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of identified targets")+xlab("Nr. of reads")
ggplot(sat[ rif == "0min" & repl == "Rep-1"],aes(N,reso,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(0,30)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Resolution")+xlab("Nr. of reads")
dev.off()


pdf(file = file.path(figsdir,"Sig70_rep1_rif20_saturation_analysis.pdf"))
ggplot(sat[ rif == "20min" & repl == "Rep-1"],aes(N,nPeaks,colour = seq))+
  geom_line(size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of candidate regions")+xlab("Nr. of reads") 
ggplot(sat[ rif == "20min" & repl == "Rep-1"],aes(N,nSites,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(380,560)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of predicted events")+xlab("Nr. of reads")
ggplot(sat[ rif == "20min" & repl == "Rep-1"],aes(N,nIden,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(55,120)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of identified targets")+xlab("Nr. of reads")
ggplot(sat[ rif == "20min" & repl == "Rep-1"],aes(N,reso,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(0,30)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Resolution")+xlab("Nr. of reads")
dev.off()




pdf(file = file.path(figsdir,"Sig70_rep2_rif0_saturation_analysis.pdf"))
ggplot(sat[ rif == "0min" & repl == "Rep-2"],aes(N,nPeaks,colour = seq))+
  geom_line(size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of candidate regions")+xlab("Nr. of reads") 
ggplot(sat[ rif == "0min" & repl == "Rep-2"],aes(N,nSites,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(380,560)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of predicted events")+xlab("Nr. of reads")
ggplot(sat[ rif == "0min" & repl == "Rep-2"],aes(N,nIden,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(55,120)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of identified targets")+xlab("Nr. of reads")
ggplot(sat[ rif == "0min" & repl == "Rep-2"],aes(N,reso,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(0,30)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Resolution")+xlab("Nr. of reads")
dev.off()


pdf(file = file.path(figsdir,"Sig70_rep2_rif20_saturation_analysis.pdf"))
ggplot(sat[ rif == "20min" & repl == "Rep-2"],aes(N,nPeaks,colour = seq))+
  geom_line(size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of candidate regions")+xlab("Nr. of reads") 
ggplot(sat[ rif == "20min" & repl == "Rep-2"],aes(N,nSites,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(380,560)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of predicted events")+xlab("Nr. of reads")
ggplot(sat[ rif == "20min" & repl == "Rep-2"],aes(N,nIden,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(55,120)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Number of identified targets")+xlab("Nr. of reads")
ggplot(sat[ rif == "20min" & repl == "Rep-2"],aes(N,reso,colour = seq))+
  stat_smooth(method = "lm",se = FALSE,formula = y ~ ns(x,3),size = 1.5)+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  theme_bw()+
  ylim(0,30)+
  theme(legend.position ="top",axis.text = element_text(size = 20),axis.title = element_text(size = 22))+
  scale_x_continuous(labels = scales::comma,limits = c(1e5,9e5))+
  ylab("Resolution")+xlab("Nr. of reads")
dev.off()

