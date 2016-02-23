rm(list = ls())

library(mosaics)
library(ChIPUtils)
library(data.table)
library(parallel)
library(devtools)
library(dpeak)


frag_len <- 150
bin_size <- 150
fdr <- .05
Gstar <- 5
mc <- detectCores()

in_dir1 <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/rif_treatment"
in_dir2 <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic"
out_dir <- "/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo"
peak_dir <- file.path(out_dir,"peaks",paste0("FDR",fdr*100))

files1 <- c("edsn1311_Sig70.sort.bam","edsn1314_Sig70.sort.bam",
           "edsn1317_Sig70.sort.bam","edsn1320_Sig70.sort.bam")
files2 <- c( "edsn931_Sig70.sort.bam","edsn933_Sig70.sort.bam")

peakfiles <- gsub(".sort.bam","_peaks.txt",c(files1,files2))

get_binding_events <- function(peakfile,readfile,Gstar,pet,frag_len,mc){

  dpeak_data <- dpeakRead(peakfile = peakfile,
    readfile = readfile,PET = pet,fragLen = frag_len,
    fileFormat = "bam",nCore = mc)
  fit <- dpeakFit(dpeak_data,nCore = mc,maxComp = Gstar)
  tt <- tempfile(pattern = "bed")
  export(fit,type = "bed",filename = tt)
  out <- data.table(read.table(tt,skip = 1))
  setnames(out,names(out),c("chrID","siteStart","siteEnd","siteName","siteStrength"))
  return(out)

}

sites <- mapply(get_binding_events,peakfile = file.path(peak_dir,peakfiles),
  readfile = c(file.path(in_dir1,files1),file.path(in_dir2,files2)),
  MoreArgs = list(Gstar = Gstar,pet = FALSE,frag_len = frag_len,
  mc = mc),SIMPLIFY = FALSE)

FF <- paste0("FDR",fdr*100)
out_dir <- file.path(out_dir,"binding_sites",FF)

dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

sitefiles <- c(gsub(".sort.bam",paste0("_sites_G",Gstar,".txt")  , files1),
               gsub(".sort.bam",paste0("_sites_G",Gstar,".txt")  , files2))
               

a = mapply(write.table,sites,file.path(out_dir,sitefiles),
  MoreArgs = list(quote = FALSE,sep = "\t", row.names = FALSE,col.names = TRUE))


