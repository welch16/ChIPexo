
rm(list = ls())
library(ChIPQC)

## Quality controls for ChIP exo data

dr ="/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo3/rawdata"
folder = c("ChIPexo","ChIPseq_PET","ChIPseq_SET")

## Set bam files on list
files = lapply(folder,function(x,dr){
  ff = list.files(file.path(dr,x))
  return(ff[!grepl("bai",ff) & !grepl("sam",ff)])},dr)
names(files) = folder

qcsamples = list()

qcsamples[["ChIPexo"]] = lapply(files[["ChIPexo"]],function(x,dr,folders)
  ChIPQCsample(file.path(dr,folders[1],x),annotation = NULL),dr,folder)
names(qcsamples[["ChIPexo"]]) = files[["ChIPexo"]]

qcsamples[["ChIPseq_PET"]] = lapply(files[["ChIPseq_PET"]],function(x,dr,folders)
  ChIPQCsample(file.path(dr,folders[2],x),annotation = NULL),dr,folder)
names(qcsamples[["ChIPseq_PET"]]) = files[["ChIPseq_PET"]]

qcsamples[["ChIPseq_SET"]] = lapply(files[["ChIPseq_SET"]],function(x,dr,folders)
  ChIPQCsample(file.path(dr,folders[3],x),annotation = NULL),dr,folder)
names(qcsamples[["ChIPseq_SET"]]) = files[["ChIPseq_SET"]]


ma <- function(x,n) filter(x,rep(1/n,n),sides = 2)

cc_simple <- function(chipqc,main)
{
  rf = QCmetrics(chipqc)[5:6]
  readl = rf[1]
  fragl = rf[2]
  cc = plotCC(chipqc)$data
  plot(cc,main = main,type = "l")
  abline(v=readl,col = "blue")
  abline(v= fragl,col = "red")
}

pdf(file = "../Figs/ChIP-Exo-CrossCorrelation.pdf")
lapply(1:length(qcsamples[["ChIPexo"]]),function(i,samp)
  cc_simple(samp[[i]],names(samp)[i]),qcsamples[["ChIPexo"]])
dev.off()

pdf(file = "../Figs/ChIP-Seq-PET-CrossCorrelation.pdf")
lapply(1:length(qcsamples[["ChIPseq_PET"]]),function(i,samp)
  cc_simple(samp[[i]],names(samp)[i]),qcsamples[["ChIPseq_PET"]])
dev.off()

pdf(file = "../Figs/ChIP-Seq-SET-CrossCorrelation.pdf")
lapply(1:length(qcsamples[["ChIPseq_SET"]]),function(i,samp)
  cc_simple(samp[[i]],names(samp)[i]),qcsamples[["ChIPseq_SET"]])
dev.off()

save(list = c("qcsamples","files","folder","dr"),file = "ChIPQC_out.RData")
