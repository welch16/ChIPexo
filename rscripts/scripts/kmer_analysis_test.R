#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    kmers_analysis_test.R 

  Arguments:


   -- exofile

      Name of the ChIP-exo file use to build the ExoData object

   -- genome


   -- mc.cores

      Number of workers used when parallel

      
   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 2)


library(ChIPexoQual)

exofile = "/p/keles/ChIPexo/volume4/carroll_data/human/bamfiles/ERR336933.sort.bam"
mc.cores = 12

exo = ExoData(exofile,mc.cores = 20,save_reads = TRUE,verbose = TRUE)

reads = exo@reads

## split by FSR

exo1 = subset(exo, FSR >= .25 & FSR <= .75)
exo2 = subset(exo, FSR < .25  | FSR > .75)

reads1 = subsetByOverlaps(reads,exo1)
reads2 = subsetByOverlaps(reads,exo2)

reads1 = split(reads1,seqnames(reads1))
reads2 = split(reads2,seqnames(reads2))

library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

seqs1 = getSeq(Hsapiens,reads1)
seqs2 = getSeq(Hsapiens,reads2)

library(magrittr)

readLength1 = sapply(width(reads1),min) %>% min
readLength2 = sapply(width(reads2),min) %>% min

readLength = min(readLength1,readLength2)

k = 4

kmerCountVec <- function(chrSeqs,initPos,k)
{
  stopifnot(initPos > 0)
  seqs = subseq(chrSeqs,start = initPos,end = initPos + k - 1)
  seqs = factor(as.character(seqs),levels = tcR::generate.kmers(k))  
  seqs %>% table
}

register(MulticoreParam(workers = mc.cores))

library(dplyr)
library(BiocParallel)

chrMats1 = lapply(seqs1[1:3],function(x,k,rl){
  positions = seq_len(rl - k + 1)
  out = bplapply(positions,
    function(y,k)kmerCountVec(x,y,k),k)
  do.call(rbind,out)  
},k,readLength) ## for testing purpuses

chrMats2 = lapply(seqs2[1:3],function(x,k,rl){
  positions = seq_len(rl - k + 1)
  out = bplapply(positions,
    function(y,k)kmerCountVec(x,y,k),k)
  do.call(rbind,out)  
},k,readLength) ## for testing purpuses


chr1 = names(chrMats1)
chr2 = names(chrMats2)

mat1 = Reduce("+",chrMats1[!chr1 %in% c("chrX","chrY","chrM")])
mat1 = mat1 / unique(rowSums(mat1))

mat2 = Reduce("+",chrMats2[!chr2 %in% c("chrX","chrY","chrM")])
mat2 = mat2 / unique(rowSums(mat2))


weights1 = colMeans(mat1[-1,]) / mat1[1,]
weights2 = colMeans(mat1[-1,]) / mat2[1,]



ww1 = data_frame(kmer = names(weights1),weight = weights1,Imbalance = "no") 
ww2 = data_frame(kmer = names(weights2),weight = weights2,Imbalance = "yes") 

ww = rbind(ww1,ww2)


library(ggplot2)

theme_set(theme_bw())

pdf(width =  12,height = 4)
ww %>% ggplot(aes(kmer,weight)) + geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90,size = 3),
        legend.position = "none")+facet_grid(  Imbalance ~ . )+
  geom_abline(slope = 0,intercept = 1,linetype = 2,colour = "red")
dev.off()
