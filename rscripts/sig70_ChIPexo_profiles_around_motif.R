
rm(list = ls())

library(GenomicRanges)
library(GenomicAlignments)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(parallel)
library(data.table)
library(dplyr)

indir <- "/p/keles/ChIPexo/volume6/K12/meme/exo_fimo"
sm <- 1
window_length <- 30
mc <- detectCores()

files <- list.files(indir,recursive = TRUE,
                    full.names = TRUE,pattern = "fimo.txt")


fimo <- lapply(files,fread)
fimo <- lapply(fimo,function(x){
  setnames(x,names(x),
           c("pattern","name","start","stop","strand",
             "score","p.val","q.val","sequence"))
  x})
names(fimo) <- basename(dirname(files))

read_dir <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo"
read_files <- list.files(read_dir,recursive = TRUE,
                         full.names = TRUE,pattern = "sort.bam")
read_files <- read_files[grep("bai",read_files,invert = TRUE)]

reads <- mclapply(read_files,readGAlignments,param = NULL,mc.cores = 8)
reads <- mclapply(reads,as,"GRanges",mc.cores = 8)
names(reads) <- gsub(".sort.bam","",basename(read_files))
depths <- sapply(reads,length)


fwd <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "+")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 8)

bwd <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "-")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 8)

rm(reads)

fimo <- lapply(fimo,function(x){
  nms <- sapply(strsplit(x[,(name)],":"),function(x)x[2])
  st <- as.numeric(sapply(strsplit(nms,"_",fixed = TRUE),function(y)y[1]))
  x[, mid := ifelse(strand == "+",st + start,st + stop)]
  return(x)  
  })


profile_DT <- function(fimo,fwd,bwd,depth = 1e9)
{
  region <- fimo[,GRanges(seqnames = "U00096",
                          ranges = IRanges(start = mid,width = 1),strand)]
  region <- resize(region,width = 2 * window_length + 1, fix = "center")
  coord <- -window_length:window_length
  fwd <- fwd[region]
  bwd <- bwd[region]

  cover <- mapply(function(sgn,f,r){
    if(sgn == "+"){
      out <- f
    }else{
      out <- r
    }
    return(out) },as.character(strand(region)),fwd,bwd,SIMPLIFY = FALSE)
  
  fimo[,name2 := paste(name,1:.N,sep = "-")]
  dt_list <- lapply(cover,function(z)data.table(coord,counts = as.vector(z)))
  dt_list <- lapply(dt_list,function(z)z[,counts := 1e9 * counts / depth])
  dt_list <- mapply(function(x,y)x[,name2 := y],dt_list,fimo[,(name2)],SIMPLIFY = FALSE)
  DT <- do.call(rbind,dt_list)
  return(merge(DT,fimo,by = "name2"))
}


edsn <- gsub("_Sig70","",names(fwd))

dt_list <- mclapply(edsn,function(x){
  message(x)
  out <- lapply(fimo[grep(x,names(fimo))],
                profile_DT,fwd[[grep(x,names(fwd))]],bwd[[grep(x,names(bwd))]],depths[grep(x,names(depths))])
  return(out)},mc.cores = 8)
#names(dt_list) <- edsn

dt_list <- unlist(dt_list,recursive = FALSE)

profile_curves <- function(prof,name,minScore = NULL,maxPval = NULL,maxQval = NULL)
{
  if(!is.null(minScore)){
    prof <- filter(prof,score >= minScore)
  }
  if(!is.null(maxPval)){
    prof <- filter(prof,p.val <= maxPval)

  }
  if(!is.null(maxQval)){
    prof <- filter(prof,q.val <= maxQval)
  }

  profiles <- prof %>% group_by(coord,strand)
  profiles <- summarize(profiles,
                        aveCounts = mean(counts),
                        trimCounts = mean(counts,trim = .1),
                        medCounts = median(counts)) %>% mutate(name = name)

  return(profiles)
      
}

profiles <- mapply(profile_curves,dt_list,names(dt_list),SIMPLIFY = FALSE)

figs_dir <- "figs/supplement"

pdf(file = file.path(figs_dir,"sig70_mean_profiles.pdf"))
plots <- lapply(profiles,function(x){
  nm <- as.character(unique(select(x,name)))
  p <- ggplot(x,aes(coord,aveCounts,colour = strand))+geom_line()+
    scale_color_brewer(palette = "Set1")+theme_bw()+
    theme(legend.position = "top")+ggtitle(nm)+
    xlab("Position around motif")+ylab("Average counts")+
    geom_vline(xintercept = 0 ,linetype = 2)  
  return(p)})
u <- lapply(plots,print)
dev.off()
 
pdf(file = file.path(figs_dir,"sig70_trim_mean_profiles.pdf"))
plots <- lapply(profiles,function(x){
  nm <- as.character(unique(select(x,name)))
  p <- ggplot(x,aes(coord,trimCounts,colour = strand))+geom_line()+
    scale_color_brewer(palette = "Set1")+theme_bw()+
    theme(legend.position = "top")+ggtitle(nm)+
    xlab("Position around motif")+ylab(".1 Trimmed mean counts")+
    geom_vline(xintercept = 0 ,linetype = 2)
  return(p)})
u <- lapply(plots,print)
dev.off()

profiles <- mapply(profile_curves,dt_list,names(dt_list),MoreArgs = list(maxPval = 1e-5),SIMPLIFY = FALSE)
    
pdf(file = file.path(figs_dir,"sig70_mean_profiles_minPval.pdf"))
plots <- lapply(profiles,function(x){
  nm <- as.character(unique(select(x,name)))
  p <- ggplot(x,aes(coord,aveCounts,colour = strand))+geom_line()+
    scale_color_brewer(palette = "Set1")+theme_bw()+
    theme(legend.position = "top")+ggtitle(nm)+
    xlab("Position around motif")+ylab("Average counts")+
    geom_vline(xintercept = 0 ,linetype = 2)  
  return(p)})
u <- lapply(plots,print)
dev.off()
 
pdf(file = file.path(figs_dir,"sig70_trim_mean_profiles_minPval.pdf"))
plots <- lapply(profiles,function(x){
  nm <- as.character(unique(select(x,name)))
  p <- ggplot(x,aes(coord,trimCounts,colour = strand))+geom_line()+
    scale_color_brewer(palette = "Set1")+theme_bw()+
    theme(legend.position = "top")+ggtitle(nm)+
    xlab("Position around motif")+ylab(".1 Trimmed mean counts")+
    geom_vline(xintercept = 0 ,linetype = 2)
  return(p)})
u <- lapply(plots,print)
dev.off()


profiles <- mapply(profile_curves,dt_list,names(dt_list),MoreArgs = list(maxQval = 2.5e-2),SIMPLIFY = FALSE)
    
pdf(file = file.path(figs_dir,"sig70_mean_profiles_minQval.pdf"))
plots <- lapply(profiles,function(x){
  nm <- as.character(unique(select(x,name)))
  p <- ggplot(x,aes(coord,aveCounts,colour = strand))+geom_line()+
    scale_color_brewer(palette = "Set1")+theme_bw()+
    theme(legend.position = "top")+ggtitle(nm)+
    xlab("Position around motif")+ylab("Average counts")+
    geom_vline(xintercept = 0 ,linetype = 2)  
  return(p)})
u <- lapply(plots,print)
dev.off()
 
pdf(file = file.path(figs_dir,"sig70_trim_mean_profiles_minQval.pdf"))
plots <- lapply(profiles,function(x){
  nm <- as.character(unique(select(x,name)))
  p <- ggplot(x,aes(coord,trimCounts,colour = strand))+geom_line()+
    scale_color_brewer(palette = "Set1")+theme_bw()+
    theme(legend.position = "top")+ggtitle(nm)+
    xlab("Position around motif")+ylab(".1 Trimmed mean counts")+
    geom_vline(xintercept = 0 ,linetype = 2)


  return(p)})
u <- lapply(plots,print)
dev.off()

fimo <- lapply(fimo,function(x)x[,name2 := paste0(name,":",.N)])

sequences <- lapply(fimo,function(x)x[,.(name2,score,p.val,q.val,sequence)])

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
},sequences,names(sequences),SIMPLIFY = FALSE)

r <- brewer.pal(name = "Set1",n = 7)

pdf(file = file.path(figs_dir,"ChIP-exo_matched_motif_sequence.pdf"),
    width = 3,height = 6 )
mapply(function(x,y){
ggplot(x,aes(position,name2,fill = fill))+
  geom_tile()+
  scale_fill_manual(values = r[c(3,2,6,1)],name = "")+
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90,size = 4),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(breaks = 1:x[,max(position)])+ggtitle(y)
},seqDT,names(seqDT),SIMPLIFY = FALSE)
dev.off()



