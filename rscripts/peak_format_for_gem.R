
rm(list = ls())

library(data.table)

base_dir <- "/p/keles/ChIPexo/volume6/resolution"

files <- list.files(base_dir,recursive = TRUE)
files <- files[grep("exo",files)]
files <- files[grep("peaks",files)]

gem_dir <- "inst/gem_analysis"

outfiles <- strsplit(files,"/",fixed = TRUE)
outfiles <- sapply(outfiles,function(x){
  out <- gsub("_peaks.txt",paste0("_",x[3],"_gem.txt"),x[4])
  return(out)})

peaks <- lapply(file.path(base_dir,files),read.table)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x)x[,paste0(V1,":",V2,"-",V3)])

a <- mapply(write.table,peaks,file.path(gem_dir,outfiles),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))            
