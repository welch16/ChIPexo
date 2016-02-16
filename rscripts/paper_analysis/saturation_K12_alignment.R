
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
afiles <- afiles[grep("bed",afiles)]

annots <- lapply(file.path(annot_dir,afiles),read.table)
annots[[1]] <- IRanges(start = annots[[1]][,2],width = 1)
annots[[2]] <- IRanges(end = annots[[2]][,3],width = 1)
annots <- reduce(c(annots[[1]],annots[[2]]))
width(annots) <- 1
annots <- sort(annots)

### test example
edsn <- data.table(exo = c(1311,1314,1317,1320),seq = c(1396,1398,1400,1402))
seeds <- list.files(file.path(indir,"ChIPexo"))
seeds <- gsub("seed","",seeds)

samps <- seq(100,900,by = 100)
ext <- 25
max_dist <- 400

annots <- annots[-1]
center <- mid(annots)
annots <- IRanges(center - ext ,center + ext)


simple_saturation_analysis <- function(samp,edsn_idx,seed,edsn,all_files,annots)
{
  
  exo_files <- all_files[grep("ChIPexo",all_files)]
  pet_files <- all_files[grep("ChIPseq_PET",all_files)]
  set_files <- all_files[grep("ChIPseq_SET",all_files)]

  samptxt <- paste0("samp",samp)

  exo_files <- exo_files[grepl(samptxt,exo_files) & grepl(edsn[edsn_idx,(exo)],exo_files) & grepl(seed,exo_files)]
  pet_files <- pet_files[grepl(samptxt,pet_files) & grepl(edsn[edsn_idx,(seq)],pet_files) & grepl(seed,pet_files)]
  set_files <- set_files[grepl(samptxt,set_files) & grepl(edsn[edsn_idx,(seq)],set_files) & grepl(seed,set_files)] 
  
  peaks <- list()
  peaks[["exo"]] <- data.table(read.table(file.path(indir,exo_files[grep("peak",exo_files)])))
  peaks[["pet"]] <- data.table(read.table(file.path(indir,pet_files[grep("peak",pet_files)])))
  peaks[["set"]] <- data.table(read.table(file.path(indir,set_files[grep("peak",set_files)])))
  
  if(nrow(peaks[["exo"]]) > 1){
    peaks[["exo"]] <- peaks[["exo"]][ V8 > quantile(V8,prob = .4)]

    sites <- list()
    sites[["exo"]] <- data.table(read.table(file.path(indir,exo_files[grep("site",exo_files)]),header = TRUE))
    sites[["pet"]] <- data.table(read.table(file.path(indir,pet_files[grep("site",pet_files)]),header = TRUE))
    sites[["set"]] <- data.table(read.table(file.path(indir,set_files[grep("site",set_files)]),header = TRUE))

    peak_ranges <- lapply(peaks,function(x)x[,IRanges(start = V2,end = V3)])
    site_ranges <- lapply(sites,function(x)x[,IRanges(start = start,end = end)])

    peak_ranges[["pet"]] <- subsetByOverlaps(peak_ranges[["pet"]],peak_ranges[["exo"]])
    peak_ranges[["set"]] <- subsetByOverlaps(peak_ranges[["set"]],peak_ranges[["exo"]])

    site_ranges <- mapply(subsetByOverlaps,site_ranges,peak_ranges,SIMPLIFY = FALSE)
  
    nPeak <- sapply(peak_ranges,length)
    nSite <- sapply(site_ranges,length)
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
  }else if(nrow(peaks[["exo"]]) == 1){
    out <- NA
  }
  return(out)
}

param <- data.table(expand.grid(k = 1:4,samp = samps,seed = seeds))
saturation <- mclapply(1:nrow(param),
  function(i){
    sat <- simple_saturation_analysis(param[i,(samp)],param[i,(k)],param[i,(seed)],
         edsn,all_files,annots)
    if(!is.na(sat) && nrow(sat) > 1 ){
      return(sat)      
    }else{
      return(NULL)
    }},mc.cores = detectCores())
    

saturation <- do.call(rbind,saturation)

sat <- saturation[,{
  nPeak = mean(nPeak)
  nSite = mean(nSite)
  nIden = mean(nIden)
  resol = mean(reso)
  list(nPeaks = nPeak,nSites = nSite,nIden = nIden,reso = resol)},by = .(seq,N,condition)]

sat[, seq := plyr::mapvalues(seq, from = c("exo","pet","set"),
        to = c("ChIP-exo","PE ChIP-Seq","SE ChIP-Seq"))]

sat[,rif := ifelse(condition %in% c(1,3),"0min","20min")]
sat[,repl := ifelse(condition %in% 1:2,"Rep-1","Rep-2")]

pdf(file = file.path(figsdir,"saturation_plots_panel.pdf"))
ggplot(sat,aes(N,nPeaks,colour = seq))+geom_line()+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  facet_grid( repl ~ rif)+theme_bw()+
  theme(legend.position ="top")+
  scale_x_continuous(labels = scales::comma,limits = c(2e5,9e5))+
  ylab("Number of candidate regions")+xlab("Nr. of reads")
ggplot(sat,aes(N,nSites,colour = seq))+geom_line()+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  facet_grid( repl ~ rif)+theme_bw()+
  theme(legend.position ="top")+
  scale_x_continuous(labels = scales::comma,limits = c(2e5,9e5))+
  ylab("Number of predicted events")+xlab("Nr. of reads")
ggplot(sat,aes(N,nIden,colour = seq))+geom_line()+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  facet_grid( repl ~ rif)+theme_bw()+
  theme(legend.position ="top")+
  scale_x_continuous(labels = scales::comma,limits = c(2e5,9e5))+
  ylab("Number of identified targets")+xlab("Nr. of reads")
ggplot(sat,aes(N,reso,colour = seq))+geom_line()+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Protocol"))+
  facet_grid( repl ~ rif)+theme_bw()+
  theme(legend.position ="top")+
  scale_x_continuous(labels = scales::comma,limits = c(2e5,9e5))+
  ylab("Resolution")+xlab("Nr. of reads")
dev.off()


