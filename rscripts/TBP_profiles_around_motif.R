
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)

motif <- "TBP"

load_all("~/Desktop/Docs/Code/ChIPexoQual")

data_dir1 <- "/p/keles/ChIPexo/volume4/venters_data/sortbam"
files1 <- list.files(data_dir1,
       pattern = "sort.bam",full.names = TRUE,include.dirs = TRUE)

files1 <- files1[grep("bai",files1,invert = TRUE)]
files1 <- files1[grep(motif,files1)]

data_dir2 <- "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam"
files2 <- list.files(data_dir2,
       pattern = "sort.bam",full.names = TRUE,include.dirs = TRUE)

files2 <- files2[grep("bai",files2,invert = TRUE)]
files2 <- files2[grep(motif,files2)]

fimo_dir <- "/p/keles/ChIPexo/volume4/tbp_analysis/fimo"
fimo_files <- list.files(fimo_dir,pattern = "txt",recursive = TRUE,
                         include.dirs = TRUE,full.names = TRUE)
fimo <- mclapply(fimo_files,read.table,mc.cores = 12)
fimo <- lapply(fimo,data.table)
names(fimo) <- basename(dirname(fimo_files))

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

fimo <- lapply(fimo,function(x)x[order(seqnames,start,end)])

window_length <- 25
sm <- 1

files <- c(files1,files2)

reads <- mclapply(files,readGAlignments,param = NULL,mc.cores = 5)
reads <- mclapply(reads,as,"GRanges",mc.cores = 5)

fwd_cover <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "+")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 5)

bwd_cover <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "-")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 5)

qv <- 1
fimo2 <- fimo[grep("TBP1",names(fimo))]

fwd_regions <- lapply(fimo2,function(x,window_length,qv){
  y <- copy(x)[strand == "+"]
  y <- y[qval <= qv]
  out <- y[,GRanges(seqnames = seqnames,
                    ranges = IRanges(start = start + motifStart,
                      width = 1))]
  out <- resize(out,width =2*window_length + 1,fix = "center")
  return(out)
},window_length,qv)

bwd_regions <- lapply(fimo2,function(x,window_length,qv){
  y <- copy(x)[strand == "-"]
  y <- y[qval <= qv]
  out <- y[,GRanges(seqnames = seqnames,
                    ranges = IRanges(start = start + motifEnd,
                      width = 1))]
  out <- resize(out,width =2*window_length,fix = "center")
  return(out)
},window_length,qv)


fwd_all <- mapply(function(cover,region,wl){
  mat <- cover[region]
  mat <- lapply(mat,as.vector)
  nms <- paste0(as.character(seqnames(region)),":",
                start(region),"-",end(region))
  DT <- mcmapply(function(x,nm,wl){
    data.table(coord = -wl : wl, counts = x , name = nm)
  },mat,nms,MoreArgs = list(wl),SIMPLIFY = FALSE,mc.cores = 10)
  return(do.call(rbind,DT))
},fwd_cover,fwd_regions,MoreArgs = list(window_length),SIMPLIFY = FALSE)

bwd_all <- mapply(function(cover,region,wl){
  mat <- cover[region]
  mat <- lapply(mat,as.vector)
  nms <- paste0(as.character(seqnames(region)),":",
                start(region),"-",end(region))
  DT <- mcmapply(function(x,nm,wl){
    data.table(coord = -wl : wl, counts = x , name = nm)
  },mat,nms,MoreArgs = list(wl),SIMPLIFY = FALSE,mc.cores = 10)
  return(do.call(rbind,DT))
},bwd_cover,bwd_regions,MoreArgs = list(window_length),SIMPLIFY = FALSE)


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
                           c("Rep-1","Rep-2","Rep-3","Rep-1","Rep-2"),
                           c(rep("ChIP-exo",3),rep("ChIP-nexus",2)),
                           SIMPLIFY = FALSE))

library(RColorBrewer)

r <- brewer.pal(name = "Set1",n = 3)[1:2]
figs_dir <- "figs/for_paper"

pdf(file = file.path(figs_dir,"TBP1_profiles_around_motif.pdf"),
    width= 8,height = 7)
ggplot(DT,aes(coord,V1,colour = strand))+
  geom_line()+
    scale_color_manual(values = rev(r))+theme_bw()+
    theme(legend.position = "top")+
  xlab("Position around motif start")+ylab("Average counts")+
  facet_grid( seq + rep  ~ . ,scales = "free")
dev.off()


sequences <- lapply(fimo2,function(x)x[,.(sequenceID,pval,qval,sequence)])

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
},sequences,c("ChIP-exo_rep1","ChIP-exo_rep2","ChIP-nexus_rep1",
              "ChIP-nexus_rep2","ChIP-nexus_rep3"),SIMPLIFY = FALSE)

r <- brewer.pal(name = "Set1",n = 7)

pdf(file = file.path(figs_dir,"TBP1_matched_motif_sequence.pdf"),
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

all_fimo <- do.call(rbind,mapply(function(x,y)x[,rep := y],
  fimo,c("ChIP-exo_rep1","ChIP-exo_rep2","ChIP-nexus_rep1",
              "ChIP-nexus_rep2","ChIP-nexus_rep3"),SIMPLIFY = FALSE))                    

pdf(file = file.path(figs_dir,"TBP1_fimo_ECDF_score.pdf"))
ggplot(all_fimo,aes(score,colour = rep))+stat_ecdf(geom = "line")+
  theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = "Replicate")+
  ylab("Empirical CDF")
ggplot(all_fimo,aes(-log10(pval),colour = rep))+stat_ecdf(geom = "line")+
  theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = "Replicate")+
  ylab("Empirical CDF")
ggplot(all_fimo,aes(-log10(qval),colour = rep))+stat_ecdf(geom = "line")+
  theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = "Replicate")+
  ylab("Empirical CDF")
dev.off()
