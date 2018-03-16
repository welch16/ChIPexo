#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

## defines option parse and help caption for script
optList = list(
  make_option("--bamfile" , action = "store_true",type = "character",
              help = "File in bam format"),
  make_option("--timefile", action = "store_true",type = "character",
              default = tempfile(),
              help = "File in tsv format where the results are going to be saved"),
  make_option("--seq_depth", action = "store_true",type = "numeric",
              default = Inf,
              help = "Sequencing depth used to sample reads without replacement"),
  make_option("--mc.cores",actio = "store_true",type = "numeric",default = 8,
              help = "Number of cores used for parallel processing")
 )


opt = parse_args(OptionParser(option_list = optList))



suppressMessages(library(GenomicAlignments))
suppressMessages(library(ChIPexoQual))
suppressMessages(library(BiocParallel))
suppressMessages(library(parallel))
suppressMessages(library(tidyverse))

options(mc.cores = opt$mc.cores)

## opt$bamfile = "/p/keles/ChIPexo/volume4/zeitlinger_data/bam/sortbam/ChIPnexus_S2_Max_rep1.sort.bam"
## opt$seq_depth = 10e6

times = list()
times[[1]] = proc.time()
reads = readGAlignments(opt$bamfile,param = NULL)

times[[2]] = proc.time()

if(!is.infinite(opt$seq_depth)){

    reads = sample(reads,opt$seq_depth)
    times[[3]] = proc.time()
    
    
}else{

    times[[3]] = times[[2]]
    
}    

exo = ExoData(reads = reads,ntimes = 1e3,nregions = 1e3)

times[[4]] = proc.time()


out = tibble(
    task = c("Reading aligned reads",
             "Sampling reads",
             "Processing ChIPexoQual pipeline"),
    time = c( times[[2]][3] - times[[1]][3],
              times[[3]][3] - times[[2]][3],
              times[[4]][3] - times[[3]][3]
             ),
    cores  = opt$mc.cores)

readr::write_tsv(x = out,path = opt$timefile)

