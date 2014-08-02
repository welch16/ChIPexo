
library(GenomicAlignments)
library(gdata)
library(parallel)

# This script only converts the data from bam file format to RData

dr ="/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo3/rawdata"
folder = c("ChIPexo","ChIPseq_PET","ChIPseq_SET")

## Set bam files on list
files = lapply(folder,function(x,dr){
  ff = list.files(file.path(dr,x))
  return(ff[!grepl("bai",ff) & !grepl("sam",ff)])},dr)
names(files) = folder
ff = lapply(folder,function(x,dr,files)file.path(dr,x,files[[x]]),dr,files)
names(ff) = folder

exo = mclapply(ff[[1]],FUN = readGAlignmentsFromBam,param = NULL,mc.cores = 8)

set = mclapply(ff[[3]],FUN = readGAlignmentsFromBam,param = NULL,mc.cores = 8)
# This part is to distinguish the files as Paired End tags

pet_flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = NA,
                hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                isFirstMateRead = NA, isSecondMateRead = NA, isNotPrimaryRead = NA,
                isNotPassingQualityControls = NA, isDuplicate = NA)
pet_param = ScanBamParam(flag = pet_flag, simpleCigar = FALSE,
                 reverseComplement = FALSE, tag = character(0),
                 what = character(0))
pet = mclapply(ff[[2]][-c(11,12)],FUN = readGAlignmentsFromBam,param = pet_param,mc.cores = 8)

names(exo) = sub("_qc.sorted.bam","",files[[1]])
names(pet) = sub("_qc.sorted.bam","",files[[2]][-(11:12)])
names(set) = sub("_qc.sorted.bam","",files[[3]])

save(list = "exo",file = "data/chip.exo.RData")
save(list = "pet",file = "data/chip.seq.pet.RData")
save(list = "set",file = "data/chip.seq.set.RData")

tab = list()
tab[[1]] = read.xls("Alignment/ChIP-Exo-summary.xls",sheet = 1)[,1:7]
tab[[2]] = read.xls("Alignment/ChIP-Exo-summary.xls",sheet = 2)[,1:7]
tab[[3]] = read.xls("Alignment/ChIP-Exo-summary.xls",sheet = 3)[,1:7]
tab[[1]]$seq = "Exo"
tab[[2]]$seq = "PET"
tab[[3]]$seq = "SET"
nn = names(tab[[1]])
tab = do.call(rbind,tab)
tab = do.call(cbind,lapply(1:ncol(tab),function(i,tab){
  if(class(tab[,i]) == "factor"){
    x = as.character(tab[,i])
    x[x == ""] = NA
  }else{
    x = tab[,i]
  }
return(x)},tab))
tab = as.data.frame(tab)
names(tab) = nn
rownames(tab)=NULL

sample.info = tab
save(list = "sample.info",file = "data/sample.summary.RData")
