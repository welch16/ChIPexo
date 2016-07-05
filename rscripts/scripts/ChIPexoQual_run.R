#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    ChIPexoQual_run.R - Applies the ChIPexoQual pipeline to partition
       the aligned reads into regions and calculate a series of summary
       statistics.

  Arguments:


   -- bamfile

      File with bam format that contains the fragments of the ChIP-exo experiment.

   -- outfile

      RData file used to store the run and estimated parameters.

  -- nreads

      Number of reads to sample. If there is no value supplied, then the
      whole set of reads is going to be used.

   -- seed

      Seed to be used

   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) <= 4)
stopifnot(length(args) > 1)

u1 <- proc.time()

library(GenomicAlignments)
library(parallel)
library(devtools)
library(data.table)
library(broom)

load_all("~/Desktop/Docs/Code/ChIPexoQual")

u2 <- proc.time()

## bamfile = "/p/keles/ChIPexo/volume4/carroll_data/mouse/ERR336935.sort.bam"
## outfile = "./test.RData"
## nreads = NULL

bamfile <- args[1]
outfile <- args[2]

if(length(args) >= 3){
  nreads <- as.numeric(args[3])
}else{
  nreads <- NULL
}

if(length(args) == 4){
  set.seed(args[4])
}

mc <- detectCores()

exo <- create_exo_experiment(bamfile,nreads = nreads,
  parallel = TRUE,mc.cores = mc)
stats <- summary_stats(exo)

ext_stats <- list()
ext_stats[["stats"]] <- stats
ext_stats[["nreads"]] <- length(reads(exo))
ext_stats[["param"]] 
ext_stats[["exec.times"]] <- data.table()

u3 <- proc.time()

u = list(u2 - u1,u3 - u2)

write.table(data.table(Reduce(rbind,u)),file = file.path("data/times",
                              gsub(".RData","_time.txt",basename(outfile))),
            quote = FALSE,row.names = FALSE,col.names = TRUE,
            sep = "\t")
                                   
save(ext_stats,file = outfile)
