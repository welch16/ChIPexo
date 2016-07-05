
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)
library(ggplot2)
library(scales)
library(RColorBrewer)

## parameters for plots
window_length <- 25
sm <- 1
topM <- 5000

## fimo results

fimo_dr <- "/p/keles/ChIPexo/volume4/tbp_analysis/fimo"
files <- list.files(fimo_dr,pattern = "fimo.txt",recursive = TRUE,
                    full.names = TRUE)

fimo <- lapply(files,fread)
names(fimo) <- basename(dirname(files))

fimo <- lapply(fimo,function(x){
  setnames(x,names(x),c("motifID","sequenceID",
                        "motifStart","motifEnd","strand",
                        "score","pval","qval","sequence"))
  return(x)})

chr_from_fimo <- function(x)sapply(strsplit(as.character(x),":"),
                function(y)y[1])

start_from_fimo <- function(x){
  out <- sapply(strsplit(as.character(x),":"),function(y)y[2])
  out <- sapply(strsplit(out,"-",fixed = TRUE),function(z)z[1])
  return(as.numeric(out))
}

end_from_fimo <- function(x){
  out <- sapply(strsplit(as.character(x),":"),function(y)y[2])
  out <- sapply(strsplit(out,"-",fixed = TRUE),function(z)z[2])
  return(as.numeric(out))
}

fimo <- lapply(fimo,function(x){
  x[,seqnames := chr_from_fimo(sequenceID)]
  x[,start := start_from_fimo(sequenceID)]
  x[,end := end_from_fimo(sequenceID)]
  return(x)})

## loading and formating reads

nexusdir <- "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"
nexusfiles <- list.files(nexusdir,recursive = TRUE,full.names = TRUE,
                         pattern = "K562")
exodir <- "/p/keles/ChIPexo/volume4/venters_data/sortbam"
exofiles <- list.files(exodir,recursive = TRUE,full.names = TRUE,
                       pattern = "TBP")
readfiles <- c(nexusfiles,exofiles)
readfiles <- readfiles[grep("bai",readfiles,invert = TRUE)]
rm(nexusdir,nexusfiles,exodir,exofiles)

reads <- mclapply(readfiles,readGAlignments,param = NULL,mc.cores = 5)
reads <- mclapply(reads,as,"GRanges",mc.cores =5)

## create coverage for both strands

fwd_cover <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "+")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 5)

bwd_cover <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "-")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 5)

## select motifs with topM scores
topM <- 8000
regions <- lapply(fimo,function(x,topM){
  x <- copy(x)[,summit := ifelse(strand == "+",start + motifStart,
                  start + motifEnd)]
  out <- x[,GRanges(seqnames = seqnames ,
                    ranges = IRanges(start = summit,width =1),
                    strand = strand)]
  out <- resize(out,2*window_length + 1,fix = "center")
  idx <- sort(x$score,index.return = TRUE,decreasing = TRUE)$ix
  return(out[idx[1:topM]])  
},topM)

fwd_regions <- lapply(regions,function(x)
  subset(x,as.character(strand(x)) == "+"))
bwd_regions <- lapply(regions,function(x)
  subset(x,as.character(strand(x)) == "-"))
                                       
fwd_all <- mapply(function(cover,region){
  wl <- window_length
  mat <- cover[region]
  mat <- lapply(mat,as.vector)
  nms <- paste0(as.character(seqnames(region)),":",
                start(region),"-",end(region))
  DT <- mcmapply(function(x,nm,wl){
    data.table(coord = -wl : wl, counts = x , name = nm)
  },mat,nms,MoreArgs = list(wl),SIMPLIFY = FALSE,mc.cores = 10)
  return(do.call(rbind,DT))
},fwd_cover,fwd_regions,SIMPLIFY = FALSE)

bwd_all <- mapply(function(cover,region){
  wl <- window_length
  mat <- cover[region]
  mat <- lapply(mat,as.vector)
  nms <- paste0(as.character(seqnames(region)),":",
                start(region),"-",end(region))
  DT <- mcmapply(function(x,nm,wl){
    data.table(coord = -wl : wl, counts = x , name = nm)
  },mat,nms,MoreArgs = list(wl),SIMPLIFY = FALSE,mc.cores = 10)
  return(do.call(rbind,DT))
},bwd_cover,bwd_regions,SIMPLIFY = FALSE)


fwd_DT <- lapply(fwd_all,function(x)x[,mean(counts),by = coord])
bwd_DT <- lapply(bwd_all,function(x)x[,mean(counts),by = coord])

profile <- function(fwd,bwd,depth,repl,sour){
  fwd <- fwd[,strand := "+"]
  bwd <- bwd[,strand := "-"]
  dt <- rbind(fwd,bwd)
  dt[,V1 := 1e9 * V1 / depth ]
  dt[,rep :=repl]
  dt[,seq := sour]
  return(dt)
}

DT <- do.call(rbind,mapply(profile,fwd_DT,
  bwd_DT,sapply(reads,length),
  c("Rep-1","Rep-2","Rep-1","Rep-2","Rep-3"),
  c(rep("ChIP-nexus",2),rep("ChIP-exo",3)),
  SIMPLIFY = FALSE))

library(RColorBrewer)

r <- brewer.pal(name = "Set1",n = 3)[1:2]
figs_dir <- "figs/for_paper"

pdf(file = file.path(figs_dir,
  paste0("TBP_profiles_around_motif_",topM,".pdf")),
    width= 8,height = 7)
ggplot(DT,aes(coord,V1,colour = strand,linetype = rep))+
  geom_line()+
    scale_color_manual(values = rev(r))+theme_bw()+
    theme(legend.position = "top")+
  xlab("Position around motif start")+ylab("Average counts")+
  facet_grid( seq   ~ . ,scales = "free")
dev.off()

pdf(file = file.path(figs_dir,
  paste0("TBP_profiles_around_motif_",topM,"_.pdf")),
    width= 8,height = 7)
ggplot(DT,aes(coord,V1,colour = strand,linetype = rep))+
  geom_line()+
    scale_color_manual(values = rev(r))+theme_bw()+
    theme(legend.position = "top")+
  xlab("Position around motif start")+ylab("Average counts")+
  facet_grid( seq +rep  ~ . ,scales = "free")
dev.off()


sequences <- lapply(fimo,function(x,topM){
  x <- copy(x)[,.(sequenceID,pval,qval,score,sequence)]
  idx <- sort(x$score,decreasing = TRUE,index.return = TRUE)$ix
  out <- x[idx[1:topM]]
  return(out)},topM)


seqDT <- mapply(function(sequ,nm){  
  motifLength <- unique(sequ[,nchar(as.character(sequence))])
  chars <- lapply(1:motifLength,function(m)sequ[,substr(sequence,m,m)])
  dt_list <- mapply(function(m,chars,sequ){
    out <- copy(sequ)
    out[,position := m]
    out[,fill := chars]
    out[,name := nm]
  },1:motifLength,chars,MoreArgs = list(sequ),SIMPLIFY = FALSE)
  return(do.call(rbind,dt_list))
},sequences,c("ChIP-nexus_Rep1","ChIP-nexus_Rep2",
              "ChIP-exo_Rep1","ChIP-exo_Rep2","ChIP-exo_Rep3"),
                SIMPLIFY = FALSE)

r <- brewer.pal(name = "Set1",n = 7)

pdf(file = file.path(figs_dir,
    paste0("TBP_matched_motif_sequence_",topM,".pdf")),
    width = 3,height = 6 )
lapply(seqDT,function(x){
  mm <- x[,max(position)]
ggplot(x,aes(position,sequenceID,fill = fill))+
  geom_tile()+
  scale_fill_manual(values = r[c(3,2,6,1)],name = "")+
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(breaks = 1:mm)+facet_grid( . ~ name )
})
dev.off()

## all_fimo <- do.call(rbind,mapply(function(x,y)x[,rep := y],
##   fimo,c("ChIP-exo_rep1","ChIP-exo_rep2","ChIP-nexus_rep1",
##               "ChIP-nexus_rep2","ChIP-nexus_rep3"),SIMPLIFY = FALSE))                    

## pdf(file = file.path(figs_dir,"TBP1_fimo_ECDF_score.pdf"))
## ggplot(all_fimo,aes(score,colour = rep))+stat_ecdf(geom = "line")+
##   theme_bw()+theme(legend.position = "top")+
##   scale_color_brewer(palette = "Set1",name = "Replicate")+
##   ylab("Empirical CDF")
## ggplot(all_fimo,aes(-log10(pval),colour = rep))+stat_ecdf(geom = "line")+
##   theme_bw()+theme(legend.position = "top")+
##   scale_color_brewer(palette = "Set1",name = "Replicate")+
##   ylab("Empirical CDF")
## ggplot(all_fimo,aes(-log10(qval),colour = rep))+stat_ecdf(geom = "line")+
##   theme_bw()+theme(legend.position = "top")+
##   scale_color_brewer(palette = "Set1",name = "Replicate")+
##   ylab("Empirical CDF")
## dev.off()
