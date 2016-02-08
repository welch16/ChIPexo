
rm(list = ls())

library(GenomicAlignments)
library(rbamtools)

dr <- "/p/keles/ChIPexo/volume7/Landick/K12"

files <- list.files(dr,recursive = TRUE)
files <- files[grep("sort",files)]
files <- files[grep("rif_tre",files)]
files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("Input",files,invert = TRUE)]

extract_info <- function(file)
{
  reader <- bamReader(file,idx=TRUE)
  readMatrix <- bamCountAll(reader)
  return(readMatrix$nAligns)
}

exofiles <- files[grep("exo",files)]
nreads <- sapply(file.path(dr,exofiles),extract_info)

nreads2sample <- 1e3 * floor(nreads / 1e3)

files <- split(file.path(dr,files),1:length(files) %% 4)

sample_commands <- function(files,nreads,outdr,isPET,seed)
{
  outdr <- file.path(outdr,sapply(strsplit(dirname(files),"/",fixed = TRUE),function(x)x[8]))
  outfiles <- basename(files)
  outfiles <- gsub(".sort.bam",paste0(".sample_seed",seed,".bam"),outfiles)
  outfiles <- file.path(outdr,outfiles)
  nreads <- ifelse(isPET, nreads / 2 , nreads)
  return(paste("rscripts/scripts/sample_N_reads.R",files,outfiles,nreads,isPET,seed)
         )
}

seed <- 123
CMDS <- unlist(mapply(sample_commands,files,nreads2sample,MoreArgs = list(outdr = file.path(dr,"sampled"),
                                             isPET = c(FALSE,TRUE,FALSE),seed = seed ),SIMPLIFY = FALSE))
out <- mclapply(CMDS,system,mc.cores = 12)
