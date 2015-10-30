

rm(list = ls())

library(GenomicAlignments)
library(data.table)
library(knitr)


dir1 <- "/p/keles/ChIPexo/volume3/LandickData"
reads_dir <- list()
reads_dir[["exo"]] <- file.path(dir1,"ChIPexo")
reads_dir[["pet"]] <- file.path(dir1,"ChIPseq_PET")
reads_dir[["set"]] <- file.path(dir1,"ChIPseq_SET")

mc <- 16

data_dir <- "data/for_paper"

get_files <- function(dir)
{
  files <- list.files(dir)
  files <- files[grep("bam",files)]
  files <- files[grep("bai",files,invert = TRUE)]
  return(files)
}

files <- lapply(reads_dir,get_files)
files <- mapply(file.path,reads_dir,files,SIMPLIFY = FALSE)

depth <- lapply(files,function(x)mclapply(x,countBam,mc.cores = mc,mc.preschedule = FALSE))

format_depth <- function(depth)
{
  depth <- do.call(rbind,depth)
  depth <- data.table(depth)
  depth <- depth[,-(1:4),with = FALSE]
  depth[,records := prettyNum(records,big.mark = ",")]
  depth[, nucleotides := NULL]
  return(depth)

}

depth <- lapply(depth,format_depth)

save(depth,file = file.path(data_dir,"EColi_experiments_depth_tables.RData"))
