rm(list = ls())

library(mosaics)
library(ChIPUtils)
library(data.table)
library(parallel)
library(devtools)


load_all("~/Desktop/Docs/Code/dpeak")

frag_len <- 150
bin_size <- 150
fdr <- .05
Gstar <- 5
mc <- 24

in_dir <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPseq_PET/rif_treatment"
out_dir <- "/p/keles/ChIPexo/volume6/K12/downstream/ChIPseq_PET"
peak_dir <- file.path(out_dir,"peaks",paste0("FDR",fdr*100))

files <- c("edsn1396_Sig70.sort.bam","edsn1398_Sig70.sort.bam",
           "edsn1400_Sig70.sort.bam","edsn1402_Sig70.sort.bam")

peakfiles <- gsub(".sort.bam","_peaks.txt",files)

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
  readfile = file.path(in_dir,files),
  MoreArgs = list(Gstar = Gstar,pet = TRUE,frag_len = frag_len,
  mc = mc),SIMPLIFY = FALSE)

FF <- paste0("FDR",fdr*100)
out_dir <- file.path(out_dir,"binding_sites",FF)

dir.create(out_dir,showWarnings = FALSE,recursive = TRUE)

sitefiles <- gsub(".sort.bam",paste0("_sites_G",Gstar,".txt")  , files)

a = mapply(write.table,sites,file.path(out_dir,sitefiles),
  MoreArgs = list(quote = FALSE,sep = "\t", row.names = FALSE,col.names = TRUE))


