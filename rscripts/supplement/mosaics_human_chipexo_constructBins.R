
rm(list = ls())

library(mosaics)
library(GenomicAlignments)
library(data.table)
library(parallel)

bin_size <- 100
frag_len <- 150

dr <- "/p/keles/ChIPexo/volume4"

outdr <- "/p/keles/ChIPexo/volume6/imbalance/bins"
files <- list.files(dr,recursive = TRUE,full.names = TRUE,
                    pattern = ".bam")
files <- files[grep("bai",files,invert = TRUE)]


### ER on MCF-7
er_files <- files[grepl("carroll",files) & grepl("human",files)]
er_files <- er_files[grep("sort",er_files)]

### CTCF on HeLa
ctcf_files <- files[grepl("pugh",files) & !grepl("rhee",files)]
ctcf_files <- ctcf_files[2:4]

### TBP on K562
tbp_files <- files[grepl("venter",files) & grepl("sort",files)]
tbp_files <- c(tbp_files[grep("TBP",tbp_files)],
               tbp_files[grep("Input",tbp_files)])


files <- c(er_files,ctcf_files,tbp_files)

## build every bin files except CTCF_hg18 (have to build that one
## manually)
a <- mclapply(files,constructBins,
              fileFormat ="bam",outfileLoc = outdr,
              fragLen = frag_len,binSize = bin_size,mc.cores = 20)
              
ctcf_files <- files[grepl("pugh",files) & !grepl("rhee",files)]

ctcf <- ctcf_files[grepl("hg18",ctcf_files) & grepl("sort",ctcf_files)][1]

reads <- readGAlignments(ctcf,param = NULL)
reads <- as(reads,"GRanges")


library(rtracklayer)


chain <- import.chain(chain_file)
hg19_reads <- liftOver(reads,chain)
hg19_reads <- unlist(hg19_reads)

reads <- resize(hg19_reads,frag_len)

chr <- seqlevels(reads)
chr <- chr[grep("rand",chr,invert = TRUE)]
chr <- chr[grep("_",chr ,invert = TRUE,fixed = TRUE)]
chr <- chr[chr != "chrM"]

reads <- split(reads,as.character(seqnames(reads)))

construct_bins <- function(chr,reads,bin_size)
{
  rr <- reads[[chr]]
  maxP <- ( 1 + floor(max(end(rr)) / bin_size)) * bin_size
  bins <- GRanges(seqnames = chr,
                  ranges = IRanges(end = seq(bin_size,maxP,
                                     by = bin_size),
                    width = bin_size))
  dt <- data.table(chr , coord = start(bins) - 1,
                   tagCounts = countOverlaps(bins,rr))

  return(dt)
}

bins <- mclapply(chr,construct_bins,reads,bin_size,mc.cores = 24)
bins <- do.call(rbind,bins)

write.table(bins,file = paste0(ctcf,"_fragL",
                   frag_len,"_bin",bin_size,".txt"),
            row.names = FALSE,col.names = FALSE,quote = FALSE)
                   
