
rm(list = ls())

library(data.table)
library(parallel)
library(GenomicAlignments)
library(BSgenome.Ecoli.NCBI.20080805)
library(devtools)

load_all("~/Desktop/Docs/Code/dpeak")

read_dir <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo"
edsn <- as.character(c(seq(1311,1320,by =3),931,933))
maxG <- 5

read_files <- list.files(read_dir,recursive = TRUE)
idx <- lapply(edsn,grep,read_files)
read_files <- do.call(c,lapply(idx,function(x)read_files[x]))
read_files <- read_files[grep("sort",read_files)]
read_files <- read_files[grep("bai",read_files,invert = TRUE)]

fdr <- "FDR1/"
peak_dir <- "/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo"
peak_files <- list.files(peak_dir,recursive = TRUE)
idx <- lapply(edsn,grep,peak_files)
peak_files <- do.call(c,lapply(idx,function(x)peak_files[x]))
peak_files <- peak_files[grep("peak",peak_files)]
peak_files <- peak_files[grep(fdr,peak_files,fixed = TRUE)]

peaks <- lapply(file.path(peak_dir,peak_files),read.table)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x)x[,V1 := "NC_010473"])

mc <- detectCores()
motif_file <- "/p/keles/ChIPexo/volume6/K12/motif_dpeak/Sig70_motifs.RData"
motif_file_par <- "/p/keles/ChIPexo/volume6/K12/motif_dpeak/Sig70_motifs_par.RData"

load(motif_file_par) ## dp_motifs

#dp_motifs <- lapply(dp_motifs,correct_motif,"NC_010473","U00096")


dpeak_call <- function(peaksfile,readsfile,motif)
{
  dpeak <- dpeakRead(peakfile = peaksfile,readfile = readsfile,fileFormat ="bam",PET = FALSE,nCore = mc)
  fit <- dpeakFit( dpeak,objectMotif= motif,maxComp = maxG,nCore = mc)
  tt <- tempfile(fileext = "bed")
  export(fit,type = ".bed",filename = tt)
  dt <- data.table(read.table(tt,skip = 1))
  setnames(dt,names(dt),c("chrID","start","end","peakID","strength"))
  return(dt)
}

## file.path(peak_dir,peak_files)

dpeakEst <- mapply(dpeak_call,file.path(peak_dir,peak_files),
       file.path(read_dir,read_files),dp_motifs,SIMPLIFY = FALSE)

save(dpeakEst,file = "dpeak_sites_motif.RData")

out_dir <- "/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo/motif_init/FDR1"
dir.create(out_dir,recursive = TRUE)

out_files <- gsub("_peaks.txt","_sites.txt",basename(peak_files))

mapply(write.table,dpeakEst,file.path(out_dir,out_files),
       MoreArgs = list(quote = FALSE , sep = "\t",row.names = FALSE))


