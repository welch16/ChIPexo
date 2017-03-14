rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)
library(ggplot2)
library(RColorBrewer)
library(ChIPexoQual)

data_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files <- list.files(data_dir,pattern = "sort.bam",full.names = TRUE,include.dirs = TRUE)
files <- files[grep("bai",files,invert = TRUE)]

reads <- mclapply(files,readGAlignments,param = NULL,mc.cores = 3)
names(reads) <- gsub(".sort.bam","",basename(files))

reads <- mclapply(reads,as,"GRanges",mc.cores = 3 )

exo_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo"

fimo_files <- list.files(exo_dir,pattern = "txt",recursive = TRUE,
                         include.dirs = TRUE,full.names = TRUE)

fimo_files = fimo_files[grepl("FOXA1",fimo_files) & !grepl("peaks",fimo_files)]


fimo <- mclapply(fimo_files,read.table,mc.cores = 12)
fimo <- lapply(fimo,data.table)
names(fimo) <- basename(dirname(fimo_files))

##     both only_reg only_peak  rep
## 1: 12726      991       502 rep1
## 2:  3271      332       442 rep2
## 3:  4074       32      3049 rep3

# V1 motif id
# V2 sequence id
# V3 start
# V4 end
# V5 strand
# V6 score
# V7 p.value
# V8 q.value
# V9 sequence

fimo <- lapply(fimo,function(x){
  setnames(x,names(x),c("motifID","sequenceID",
                        "motifStart","motifEnd","strand",
                        "score","pval","qval","sequence"))
  return(x)})

fimo <- fimo[grep("FOXA1",names(fimo))]

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

fimo2 <- lapply(fimo,function(x)x[qval <= .05])

window_length <- 20
sm <- 1

## want to see if there is a spatial pre-disposition..
## we are going to divide the motifs by the strand and build a + / - 50 bp
## profile for each strand

fwd_cover <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "+")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 3)

bwd_cover <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "-")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 3)

fwd_regions <- lapply(fimo2,function(x,window_length){
  y <- copy(x)[strand == "+"]
  out <- y[,GRanges(seqnames = seqnames,
                    ranges = IRanges(start = start + motifStart,
                      width = 1))]
  out <- resize(out,width =2*window_length + 1,fix = "center")
  return(out)
},window_length)

bwd_regions <- lapply(fimo2,function(x,window_length){
  y <- copy(x)[strand == "-"]
  out <- y[,GRanges(seqnames = seqnames,
                    ranges = IRanges(start = start + motifEnd,
                      width = 1))]
  out <- resize(out,width =2*window_length,fix = "center")
  return(out)
},window_length)


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

profile <- function(fwd,bwd,depth,repl){
  fwd <- fwd[,strand := "+"]
  bwd <- bwd[,strand := "-"]
  dt <- rbind(fwd,bwd)
##  dt[,V1 := 1e9 * V1 / depth ]
  dt[,rep :=repl]
  return(dt)
}

DT <- do.call(rbind,mapply(profile,fwd_DT,
                           bwd_DT,sapply(reads,length),
                           c("Rep-3","Rep-1","Rep-2"),SIMPLIFY = FALSE))

library(dplyr)
library(readr)
library(tidyr)


DT = DT %>% as.tbl %>% rename(counts = V1)

write_tsv(DT,path = "data/figures/fig4/fig4_profiles_FoxA1.tsv")


## sequence plots

sequences <- lapply(fimo,function(x)x[,.(sequenceID,pval,qval,sequence)])

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
},sequences,c("Rep-3","Rep-1","Rep-2"),SIMPLIFY = FALSE)


seqDT = seqDT %>% lapply(as.tbl) %>% bind_rows

write_tsv(seqDT,"data/figures/fig4/fig4_sequences_FoxA1.tsv")

outdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"

## check that the "high quality" regions overlap with peaks
peakfiles <- list.files(outdir,
   pattern = "peaks",recursive = TRUE,full.names= TRUE)
peaks <- lapply(peakfiles,
               read.table)
names(peaks) <- basename(peakfiles)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x)x[,1:3,with = FALSE])
peaks <- lapply(peaks,function(x){
  setnames(x,names(x),c("seqnames","start","end"))
  return(ChIPUtils::dt2gr(x))})

regions <- lapply(fimo,
        function(x)ChIPUtils::dt2gr(x[,.(seqnames,start,end)]))


venn <- mapply(function(region,peak){
  ov <- findOverlaps(region,peak)
  inter_peak <- length(unique(subjectHits(ov)))
  inter_reg <- length(unique(queryHits(ov)))
  out <- data.table(both = inter_reg,
                    only_reg = length(region) - inter_reg,
                    only_peak = length(peak) - inter_peak)
  return(out)                   
},regions,peaks,SIMPLIFY = FALSE)
venn <- do.call(rbind,venn)

venn[, rep := c("rep3","rep1","rep2")]

write_tsv(venn,"data/figures/fig4/fig4_vennDiagrams.tsv")

options(mc.cores = 20)

exo = reads %>% lapply(function(x)ExoData(reads = x))

rl = reads %>% lapply(width) %>% lapply(median)

rl = rl %>% unlist %>% unique


stats = exo %>% lapply(function(x)x %>% as.data.frame %>% as.tbl)

## 1 - only regions with reads in both strands
stats = stats %>% lapply(function(x)x %>% filter(fwdReads > 0 & revReads > 0))

## 2 - wider regions
stats = stats %>% lapply(function(x)x %>% filter( width >= 3 * rl & width <= 1e3))

                                   
## 3 - deeper regions
minNpos <- 15
minDepth <- 100

stats = stats %>% lapply(function(x)x %>% filter(uniquePos >= minNpos) %>%
                                    filter(depth >= minDepth))

## remove chrM regions

stats = stats %>% lapply(function(x)x %>% filter(seqnames != "chrM"))

stats_gr = stats %>% lapply(function(x){
    GRanges(seqnames = as.character(x$seqnames),
            ranges = IRanges(
                start = x$start,end = x$end))})

fimo = fimo %>% lapply(as.tbl)

fimo_gr = fimo %>% lapply(function(x){
    GRanges(seqnames = as.character(x$seqnames),
            ranges = IRanges(
                start = x$start,end = x$end))})

names(fimo) = paste("Rep",c(3,1,2),sep = "-")
fimo = mapply(function(x,y)x %>% mutate(Replicate = y),fimo,names(fimo),SIMPLIFY = FALSE)


topK <- c(50,250,100,500,1000,2000)


dt_list = fimo %>% lapply(function(x){
    lapply(topK,function(z)top_n(x,z,score) %>%
                           mutate(K = z)) %>% bind_rows}) %>% bind_rows

write_tsv(dt_list,"data/figures/fig4/fig4_FIMO_scores_FoxA1.tsv")

all_peaks = Reduce(c,peaks) %>% reduce

exo_peaks = mclapply(exo,subsetByOverlaps,all_peaks)


## remove chrM
exo_peaks = exo_peaks %>% mclapply(function(x)x[as.character(seqnames(x)) != "chrM"])


K = 3
exo_peaks = exo_peaks %>% lapply(function(x)x[width(x) >= K * rl])

all_fimo = fimo %>% bind_rows


strand_imbalance_fimo <- function(exo,fimo,repl)
{
    exo_tbl = exo %>% as.data.frame %>% as.tbl %>%
        mutate(match = paste0(seqnames,":",start,"-",end)) %>%
        select(match , everything())
    fimo = fimo %>% rename(match = `sequenceID`) %>% mutate(match = as.character(match))

    fimo_summary = fimo %>% group_by(match) %>%
        summarize(
            nmotif = n(),
            minPval = min(`pval`),
            maxScore = max(score)
        )

    exo_tbl %>% left_join(fimo_summary,by ="match") %>%
        mutate(nmotif = ifelse(is.na(nmotif),0,nmotif),repl)

  
}


imbalance = mapply(strand_imbalance_fimo,exo_peaks,fimo,names(fimo),SIMPLIFY = FALSE) %>% bind_rows

maxMotif = max(imbalance$nmotif)

imbalance = imbalance %>%
    mutate(nmotif = forcats::fct_collapse(as.factor(nmotif),">=3" = as.character(seq(3,maxMotif))))

write_tsv(imbalance,"data/figures/fig4/fig4_imbalance_hist_FoxA1.tsv")

