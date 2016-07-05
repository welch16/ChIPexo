#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    SCC_curve_and_QC_ChIPseq_SE.R - Calculates the SCC curve for a
      given SE ChIP-seq sample.

  Arguments:

   -- bedfile

      File in bed format with the aligned reads used to calculate the curve

   -- outfile

      Name of the file, where the output is gonna be saved

   -- summaryfile

      Name of the file, where a collection of ChIP-Seq QC metrics are
      gonna be saved

   -- sizefile

      File without header and with two column indicating the
      respective chromosome and it's length. If the data was aligned to
      any of the dm3, hg19, mm9 or mm10 genomes, then it loads the file
      automatically by using dm3, hg19, etc.

   -- isPET

      Boolean variable indicating if the reads in the bamfile are paired

   -- maxShift

      Max possible shift for the curve, i.e. the output is gonna be a
      data.table with format [shift = 1:maxShift , cross.corr]

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 6)

bedfile <- args[1]
outfile <- args[2]
summaryfile <- args[3]
sizefile <- args[4]
isPET <- as.logical(args[5])
maxShift <- as.numeric(args[6])

stopifnot(file.exists(bedfile))
stopifnot(maxShift > 0)

library(parallel)
library(data.table)
library(GenomicAlignments)
library(devtools)

load_all("~/Desktop/Docs/Code/ChIPUtils")

mc <- detectCores()

create_reads_dt_set <- function(dt,reads_file){

  dt <- fread(reads_file)
  
  greads <- dt[,GRanges(seqnames = V1,
    ranges = IRanges(start = V2,end = V3),strand = V6)]
  
  gr <- gr2dt(greads)
  setkey(gr,strand)

  fwd <- gr["+",nomatch = 0]
  bwd <- gr["-",nomatch = 0]
  fwd <- split(fwd,fwd[,(seqnames)])
  bwd <- split(bwd,bwd[,(seqnames)])

  out <- new("reads",readsFile = reads_file,readsF = fwd,readsR = bwd,
             nReads = nrow(gr),isPET = FALSE)
  return(out)
}

reads <- create_reads_dt_set(dt,bedfile)


if(tolower(sizefile) %in% c("hg19","mm9","mm10","dm3")){
  sizedir <- system.file("extdata","chrom.sizes", package = "ChIPUtils")
  sizefiles <- list.files(sizedir)
  sizefile <- sizefiles[grep(sizefile,sizefiles)]
  sizefile <- file.path(sizedir,sizefile)
  rm(sizedir,sizefiles)
}

sizes <- data.table(read.table(sizefile,header = FALSE))

scc <- strand_cross_corr(reads,shift = 1:maxShift,
   chrom.sizes = sizes,parallel = FALSE)

write.table(format(scc,digits = 6),file = outfile,quote = FALSE,
   sep = "\t",row.names = FALSE,col.names = TRUE)            

strand_ratio <- function(reads){
  fwd <- length(dt2gr(do.call(rbind,readsF(reads))))
  bwd <- length(dt2gr(do.call(rbind,readsR(reads))))
  out <- fwd / (fwd + bwd)
  return(out)
}
  
read_length <- function(reads){
  fwd <- dt2gr(do.call(rbind,readsF(reads)))
  bwd <- dt2gr(do.call(rbind,readsR(reads)))
  out <- c(fwd,bwd)
  out <- table(width(out))
  out <- as.numeric(names(which.max(out)))
  return(out)
}


NSC <- function(scc)scc[,max(cross.corr) /min(cross.corr)]

RSC1 <- function(scc,read_length){
  out <- scc[,max(cross.corr)] / scc[shift == read_length, (cross.corr)]
  return(out)
}

RSC2 <- function(scc,read_length){
  mm <- scc[,min(cross.corr)]
  out <- (scc[,max(cross.corr)] - mm) /( scc[shift == read_length, (cross.corr)] - mm)
  return(out)
}


rl <- read_length(reads)
fl <- scc[which.max(cross.corr),(shift)]

summary <- data.table(depth = nreads(reads),
                      PBC = PBC(reads),
                      FSR = strand_ratio(reads),
                      read_length = rl,
                      frag_length = fl,
                      NSC =  NSC(scc),
                      RSC1 = RSC1(scc,rl),
                      RSC2 = RSC2(scc,rl))
                      
write.table(format(summary,digits = 6),file = summaryfile,quote = FALSE,
   sep = "\t",row.names = FALSE,col.names = TRUE)            
