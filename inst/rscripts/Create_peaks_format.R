
rm(list = ls())
library(data.table)

extract_transform_save <- function(filecodename,filename,basedir)
{
  df = read.table(file = file.path(basedir,filename),stringsAsFactors= FALSE,header=FALSE)
  peaks = data.table(chrID = df[,1],peakStart = df[,2],peakStop = df[,3],score = df[,5])
  savedir = file.path(basedir,filecodename,"peaks")
  if(!file.exists(savedir))dir.create(savedir)
  message(file.path(savedir,paste0(filecodename,"-peaks.RData")))
  save(peaks,file = file.path(savedir,paste0(filecodename,"-peaks.RData")))
}



# Carroll human
dir = "/p/keles/ChIPexo/volume3/Analysis/Carroll/human"
datasets = paste0("ER-rep",1:3)

## files = c("jc1199_MCF7_ER_exo_repA_CRI01_hg19_MACS2_jc899_peaks.bed",
##           "jc1200_MCF7_ER_exo_repB_CRI01_hg19_MACS2_jc899_peaks.bed",  
##           "jc1201_MCF7_ER_exo_repC_CRI01_hg19_MACS2_jc899_peaks.bed")

files = c("jc1202_MCF7_ER_ChIP-seq_repA_hg19_MACS2_jc899_peaks.bed",
          "jc1203_MCF7_ER_ChIP-seq_repB_hg19_MACS2_jc899_peaks.bed",
           "jc1204_MCF7_ER_ChIP-seq_repC_hg19_MACS2_jc899_peaks.bed")


  ## 37211 jc1199 rep 1 exo
  ## 43960 jc1200 rep 2 exo
  ## 46151 jc1201 rep 3 exo
  ## 27501 jc1202 rep 1 seq
  ## 26661 jc1203 rep 2 seq
  ## 24234 jc1204 rep 3 seq

mapply(extract_transform_save,datasets,files,MoreArgs = list(basedir = dir))


# Carroll mouse (this are not peaks but mesas as they called them)
dir = "/p/keles/ChIPexo/volume3/Analysis/Carroll/mouse"
datasets = paste0("FoxA1-rep",1:3)
files = c("do2160_FoxA1_exo_mmuBL6_0SQ40_CRI01_mm9_mesa.bed",
          "do2161_FoxA1_exo_mmuBL6_0SQ41_CRI01_mm9_mesa.bed",
          "do2435_FoxA1_exo_mmuBL6_0SQ42_CRI02_mm9_mesa.bed")

mapply(extract_transform_save,datasets,files,MoreArgs = list(basedir = dir))


# Ren
dir = "/p/keles/ChIPexo/volume3/Analysis/Ren"
datasets = c("H3k27ac",paste0("H3k4me1-rep",1:2))

files = c("K562_H3k27ac.bed",rep("K562_H3k4me3_broad.bed",2))

mapply(extract_transform_save,datasets,files,MoreArgs = list(basedir = dir))



# Landick
dir = "/p/keles/ChIPexo/volume3/Analysis/Landick"
load(file.path(dir,"peaks.RData"))

dt1 = peaks_list[1:4]
dt2 = peaks_list[5:8]

dir1 = file.path(dir,"stat-vs-exp")
dir2 = file.path(dir,"rif-treatment")


filenames1 = names(dt1)
filenames2 = names(dt2)

filecodenames1 = c("Sig70-exp-rep1",
                   "Sig70-exp-rep2",
                   "Sig70-stat-rep1",
                   "Sig70-stat-rep2")
filecodenames2 = c("Sig70-rif0-rep1",
                   "Sig70-rif20-rep1",  
                   "Sig70-rif0-rep2",
                   "Sig70-rif20-rep2")

extract_transform_save <- function(filecodename,dt,basedir)
{
  peaks = data.table(chrID = dt$chrID,peakStart = dt$peakStart,peakStop = dt$peakStop,score = dt$aveP)
  savedir = file.path(basedir,filecodename,"peaks")
  if(!file.exists(savedir))dir.create(savedir)
  message(file.path(savedir,paste0(filecodename,"-peaks.RData")))
  save(peaks,file = file.path(savedir,paste0(filecodename,"-peaks.RData")))
}

mapply(extract_transform_save,filecodenames1,dt1,MoreArgs = list(dir1))
mapply(extract_transform_save,filecodenames2,dt2,MoreArgs = list(dir2))
