
rm(list = ls())

library(data.table)
library(parallel)

frag_len <- 150
bin_size <- 100

dr <- "/p/keles/ChIPexo/volume6/imbalance/addfiles"

dirs <- list.files(dr,include.dirs = TRUE,full.names = TRUE)
names(dirs) <- basename(dirs)

files <- lapply(dirs,list.files,full.names = TRUE)
files <- lapply(files,function(x)x[grep("temp",x,invert = TRUE)])

chrs <- lapply(files,function(x){
  nms <- basename(x)
  out <- sapply(strsplit(nms,"_"),function(y)y[1])
  return(out)})

dts <- lapply(files,function(x){
  out <- mclapply(x,fread,mc.cores = 20)
  return(out)})
              
dts <- mapply(function(x,y){
  out <- mcmapply(function(z,w){
    z[,chr := w]
    return(z)},x,y,SIMPLIFY = FALSE,mc.cores = 20)
  return(do.call(rbind,out))},dts,chrs,SIMPLIFY = FALSE)

dts <- lapply(dts,function(x){
  setnames(x,names(x),c("coord","val","chr"))
  setcolorder(x,c("chr","coord","val"))
  return(x)})

outfiles <- file.path(dr,paste0(names(dts),"_fragL",frag_len,"_binSize",
                                bin_size,".txt"))

mcmapply(write.table,dts,outfiles,
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE),
         mc.cores = 3 )       




