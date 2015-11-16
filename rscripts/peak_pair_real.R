
rm(list = ls())

library(data.table)
library(GenomicRanges)
library(GenomicAlignments)
library(viridis)
library(scales)
library(ggplot2)
library(ChIPUtils)

base_dir <- "/p/keles/ChIPexo/volume3/LandickData"

exo_dir <- file.path(base_dir,"ChIPexo")
pet_dir <- file.path(base_dir,"ChIPseq_PET")
set_dir <- file.path(base_dir,"ChIPseq_PET")


## peaks
exo_peaks <- "/p/keles/ChIPexo/volume6/peak_pair/peaks/edsn931_042814_peaks.txt"
pet_peaks <- "/p/keles/ChIPexo/volume6/peak_pair/peaks/edsn790_042814_filter_peaks.txt"
set_peaks <- "/p/keles/ChIPexo/volume6/peak_pair/peaks/edsn80_042814_peaks.txt"

exo_peaks <- data.table(read.table(exo_peaks))
pet_peaks <- data.table(read.table(pet_peaks))
set_peaks <- data.table(read.table(set_peaks))

exo_gr <- exo_peaks[,1:3,with = FALSE]
pet_gr <- pet_peaks[,1:3,with = FALSE]
set_gr <- set_peaks[,1:3,with = FALSE]

setnames(exo_gr,names(exo_gr),c("seqnames","start","end"))
setnames(pet_gr,names(pet_gr),c("seqnames","start","end"))
setnames(set_gr,names(set_gr),c("seqnames","start","end"))

exo_gr <- dt2gr(exo_gr)
pet_gr <- dt2gr(pet_gr)
set_gr <- dt2gr(set_gr)

common <- subsetByOverlaps(exo_gr, subsetByOverlaps(pet_gr,set_gr))

exo_file <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo/edsn931_042814_qc.sorted.bam"
pet_file <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET/edsn790_042814_qc.filter.bam"
set_file <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET/edsn80_042814_qc.sorted.bam"
 
prm <- NULL

exo <- readGAlignments(exo_file, param = prm)
pet <- readGAlignmentPairs(pet_file,param = prm)
set <- readGAlignments(set_file,param = prm)

## fix pet
l <- left(pet)
r <- right(pet)

c1 <- start(l)
c2 <- end(r)

idx <- c1 < c2
if(any(idx)){
  s <- seqnames(pet)[idx]
  pet <- GRanges(seqnames = s,
                 ranges = IRanges(start = c1[idx],
                   end = c2[idx]),strand = as.character(strand(pet))[idx])
}

exo <- as(exo,"GRanges")
set <- as(set,"GRanges")

ma_plot_dt <- function(reads,peaks,bin_size = 150, M )
{
  n0 <- length( reads)
  message(prettyNum(n0,big.mark = ","))
  ov <- findOverlaps(reads,peaks)
  n1 <- length(ov)
  message(prettyNum(n1,big.mark = ","))
  message(prettyNum(n1 / n0 * 100) )

  bins <- create_bins(bin_size, chrom = GRanges("U00096",
    ranges = IRanges(start = 1, end=  max(end(reads)))))

  DT1 <- gr2dt(bins)
  DT2 <- copy(DT1)

  DT1[  , f := countOverlaps(bins,subset(reads,strand(reads) == "+"))]
  DT1[  , r := countOverlaps(bins,subset(reads,strand(reads) == "-"))]

  DT1 <- DT1[ f > 0 & r > 0 ]
  DT1[, label := "All bins"]

  reads_in_peaks <- reads[queryHits(ov)]

  DT2[ , f := countOverlaps(bins,subset(reads_in_peaks,strand(reads_in_peaks) == "+"))]
  DT2[ , r := countOverlaps(bins,subset(reads_in_peaks,strand(reads_in_peaks) == "-"))]

  DT2 <- DT2[ f > 0 & r > 0]
  DT2[ , label := "Bins in peaks" ]

  DT <- rbind(DT1,DT2)

  return(DT)
}

bs <- 150
M <- max(max(end(exo)), max(end(pet)),max(end(set)))
exo_DT <- ma_plot_dt(exo,bin_size = bs,common,M = M )
pet_DT <- ma_plot_dt(pet,bin_size = bs,common,M = M )
set_DT <- ma_plot_dt(set,bin_size = bs,common,M = M )

exo_DT[,what := "exo"]
pet_DT[,what := "pet"]
set_DT[,what := "set"]
DT <- rbind(exo_DT,pet_DT,set_DT)

DT[ , A := log2(f / r)]
DT[ , M := log2( f * r)]

DT[ , fsr := f / (f + r)]

ggplot(DT[what != "pet" & grepl("All",label)] , aes(fsr,colour = what)) + geom_density()+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  scale_color_brewer(palette = "Set1")+ylab("")+geom_vline( xintercept = .5 ,linetype =2 )+
  xlab("")
dev.off()







