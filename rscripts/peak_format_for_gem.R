
rm(list = ls())
library(data.table)

basedir <- "/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo"

files <- list.files(basedir,recursive = TRUE)
files <- files[grep("peaks",files)]

gem_dir <- "/p/keles/ChIPexo/volume6/K12/other_methods/gem/peaks"

outfiles <- basename(files)

sep_files <- strsplit(files,"/",fixed = TRUE)
fdr <- sapply(sep_files,function(x)x[2])

outfiles <- mapply(function(infile,fdr){
  out <- gsub("_peaks",paste0("_peaks_",fdr),infile)
  out},outfiles,fdr)
names(outfiles) <- NULL

peaks <- lapply(file.path(basedir,files),read.table)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x)x[,paste0(V1,":",V2,"-",V3)])

a <- mapply(write.table,peaks,file.path(gem_dir,outfiles),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))            
