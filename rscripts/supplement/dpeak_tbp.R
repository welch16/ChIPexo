
rm(list = ls())

library(GenomicAlignments)
library(data.table)
library(parallel)

bin_size <- 100
frag_len <- 150

dr <- "/p/keles/ChIPexo/volume4/venters_data"
files <- list.files(dr,pattern = "bam",recursive = TRUE,
                    full.names = TRUE)
files <- files[grep("bai",files,invert = TRUE)]
files <- files[grep("input",tolower(files),invert = TRUE)]
files <- files[grep("sort",files)]
files <- files[grep("TBP",files)]

devtools::load_all("~/Desktop/Docs/Code/dpeak")

load("data/TBP_peaks.RData")

lapply(peaks,nrow)

save_peaks <- function(x,topK = 1500){
  tt <- tempfile(pattern = "peaks",fileext = ".txt")

  write.table(x[order(-aveChipCount)][1:topK],
              file = tt,col.names = FALSE,
             row.names = FALSE,quote = FALSE)

  return(tt)

}

pp <- lapply(peaks,save_peaks,topK = 1500)

outdr <- "/p/keles/ChIPexo/volume6/imbalance/binding_sites"

dpeak_model <- function(read_file,peak_file,outname){
  dpeak <- dpeakRead(peakfile = peak_file,readfile = read_file,
                     fileFormat = "bam",nCore = 20)
  fit <- dpeakFit(dpeak,nCore = 20,maxComp = 5)
  export(fit,filename = outname)

}

dpeak_model(files[1],pp[[4]],
    file.path(outdr,"TBP_chipseq_sites_rep1.txt"))
dpeak_model(files[2],pp[[5]],
    file.path(outdr,"TBP_chipseq_sites_rep2.txt"))
dpeak_model(files[3],pp[[1]],
    file.path(outdr,"TBP_chipexo_sites_rep1.txt"))
dpeak_model(files[4],pp[[2]],
    file.path(outdr,"TBP_chipexo_sites_rep2.txt"))
dpeak_model(files[5],pp[[3]],
    file.path(outdr,"TBP_chipexo_sites_rep3.txt"))

