
rm(list = ls())

library(GenomicAlignments)
library(GenomicRanges)
library(data.table)
library(parallel)

fl <- 1

dr <- "/p/keles/ChIPexo/volume6/imbalance/binding_sites"
files <- list.files(dr,full.names = TRUE)

ctcf <- files[grep("ctcf",files)]
tbp <- files[grep("TBP",files)]
er <- files[grep("ER",files)]

load_sites <- function(files)
{
  sites <- lapply(files,fread)
  names(sites) <- gsub(".txt","",basename(files))
  return(sites)
}

ctcf_sites <- load_sites(ctcf)
tbp_sites <- load_sites(tbp)
er_sites <- load_sites(er) 

get_peaks <- function(sites)
{
  x <- strsplit(unique(sites$V4),":")
  sqnms <- sapply(x,function(y)y[1])
  x <- strsplit(sapply(x,function(y)y[2]),"-")
  ra <- IRanges(start = as.numeric(sapply(x,function(y)y[1])),
                end = as.numeric(sapply(x,function(y)y[2])))
  out <- GRanges(seqnames = sqnms,ranges = ra)
  sites <- sites[,GRanges(seqnames = V1,
        ranges = IRanges(start =
          mid(IRanges(start = V2,end = V3)),width = 1))]
  out$nsites <- countOverlaps(out,sites)
  return(out)
}

strand_fragments <- function(peaks,reads,side)
{
  out <- countOverlaps(peaks,
    subset(reads,strand(reads) == side))
  return(out)
}


ctcf_peaks <- lapply(ctcf_sites,get_peaks)
tbp_peaks <- lapply(tbp_sites,get_peaks)
er_peaks <- lapply(er_sites,get_peaks)

reads <- list.files("/p/keles/ChIPexo/volume4",full.names = TRUE,
                    recursive = TRUE,pattern = "bam")
reads <- reads[grep("bai",reads,invert = TRUE)]

#### CTCF

ctcf_reads_files <- reads[grepl("ctcf",tolower(reads))]
ctcf_reads_files <- ctcf_reads_files[c(1,2,4)]
ctcf_reads <- mclapply(ctcf_reads_files,
       readGAlignments,param = NULL, mc.cores = 10)
ctcf_reads <- mclapply(ctcf_reads,as,"GRanges",mc.cores = 10)
ctcf_reads <- mclapply(ctcf_reads,resize,fl,mc.cores = 10)
names(ctcf_reads) <- gsub(".bam","",basename(ctcf_reads_files))

ctcf_peaks <- mcmapply(function(peaks,reads){
  peaks$fwd <- strand_fragments(peaks,reads,"+")
  peaks$bwd <- strand_fragments(peaks,reads,"-")
  return(peaks)},ctcf_peaks,ctcf_reads[c(3,1:2)],
   mc.cores = 10,SIMPLIFY = FALSE)

save(ctcf_peaks,file = "data/CTCF_strand_imbalance.RData")

#### TBP

tbp_reads_files <- reads[grepl("tbp",tolower(reads))]
tbp_reads_files <- tbp_reads_files[grep("subsam",
          tbp_reads_files,invert = TRUE)]
tbp_reads_files <- tbp_reads_files[grepl("sort",tbp_reads_files) &
   ! grepl("zeit",tbp_reads_files)]

tbp_reads <- mclapply(tbp_reads_files,
       readGAlignments,param = NULL, mc.cores = 10)
tbp_reads <- mclapply(tbp_reads,as,"GRanges",mc.cores = 10)
tbp_reads <- mclapply(tbp_reads,resize,fl,mc.cores = 10)
names(tbp_reads) <- gsub(".bam","",basename(tbp_reads_files))

tbp_peaks <- mcmapply(function(peaks,reads){
  peaks$fwd <- strand_fragments(peaks,reads,"+")
  peaks$bwd <- strand_fragments(peaks,reads,"-")
  return(peaks)},tbp_peaks,tbp_reads[c(3:5,1:2)],
   mc.cores = 10,SIMPLIFY = FALSE)

save(tbp_peaks,file = "data/TBP_strand_imbalance.RData")

#### ER

er_reads_files <- reads[grepl("carr",tolower(reads)) &
                        grepl("hum",tolower(reads))]
er_reads_files <- er_reads_files[grepl("sort",er_reads_files)]
er_reads_files <- er_reads_files[grep("input",er_reads_files,
                                      invert = TRUE)]

er_reads <- mclapply(er_reads_files,
       readGAlignments,param = NULL, mc.cores = 10)
er_reads <- mclapply(er_reads,as,"GRanges",mc.cores = 10)
er_reads <- mclapply(er_reads,resize,fl,mc.cores = 10)
names(er_reads) <- gsub(".bam","",basename(er_reads_files))

er_peaks <- mcmapply(function(peaks,reads){
  peaks$fwd <- strand_fragments(peaks,reads,"+")
  peaks$bwd <- strand_fragments(peaks,reads,"-")
  return(peaks)},er_peaks,er_reads[c(4:6,1:3)],
   mc.cores = 10,SIMPLIFY = FALSE)

save(er_peaks,file = "data/ER_strand_imbalance.RData")

