
rm(list = ls())

library(dpeak)
library(data.table)
library(parallel)
library(BSgenome.Ecoli.NCBI.20080805)

## read_dir <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo"
edsn <- as.character(c(seq(1311,1320,by =3),931,933))
mc <- 22

## read_files <- list.files(read_dir,recursive = TRUE)
## idx <- lapply(edsn,grep,read_files)
## read_files <- do.call(c,lapply(idx,function(x)read_files[x]))
## read_files <- read_files[grep("sort",read_files)]
## read_files <- read_files[grep("bai",read_files,invert = TRUE)]

fdr <- "FDR5/"
peak_dir <- "/p/keles/ChIPexo/volume6/K12/downstream/ChIPexo"
peak_files <- list.files(peak_dir,recursive = TRUE)
idx <- lapply(edsn,grep,peak_files)
peak_files <- do.call(c,lapply(idx,function(x)peak_files[x]))
peak_files <- peak_files[grep("peak",peak_files)]
peak_files <- peak_files[grep(fdr,peak_files,fixed = TRUE)]

peaks <- lapply(file.path(peak_dir,peak_files),read.table)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x)x[,V1 := "NC_010473"])

## topM <- 500
## peaks <- lapply(peaks,function(x)x[order(-V8)][1:topM])

peak_tmp <- replicate(length(peaks),{tempfile(pattern = "peak_",fileext = ".txt")})
mapply(write.table,peaks,peak_tmp,MoreArgs = list(quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE))

memeArg <- paste("-dna -mod zoops -nmotifs 1 -minw 10 -maxw 20 -revcomp -maxsize 1000000000 -nostatus","-p",mc,sep = " ")
fimoArg <- "-max-stored-scores 100000000 -motif-pseudo 0.000001 --verbosity 1"

motif_dir <- "/p/keles/ChIPexo/volume6/K12/motif_dpeak"

dp_motifs <- lapply(peak_tmp,dpeakMotif,refGenome=Ecoli,memeArgument = memeArg,fimoArgument = fimoArg )

## motif_dir <- file.path(motif_dir,paste0("top",topM))
## dir.create(motif_dir,recursive = TRUE)

names(dp_motifs) <- basename(peak_files)

save(dp_motifs,file = file.path(motif_dir,"Sig70_motifs_FDR5.RData"))
