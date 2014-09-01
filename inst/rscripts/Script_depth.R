
rm(list = ls())
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)
library(reshape2)
library(WriteXLS)


load("data/chip.exo.RData")
load("data/chip.seq.pet.RData")
load("data/chip.seq.set.RData")
load("data/sample.summary.RData")

source("R/density_functions.R")
source("R/conditional_density_functions.R")
source("R/depth_functions.R")

mc.cores = 8
binsizes = c(200,500,750,1000)
figsdir = "inst/figs/depth"

exo = mclapply(exo,as.GRanges,mc.cores=mc.cores)
pet = mclapply(pet,as.GRanges,mc.cores=mc.cores)
set = mclapply(set,as.GRanges,mc.cores=mc.cores)

index = !grepl("93?",names(exo))
exo_sets1 = exo[index]
exo_sets2 = exo[!index]  

# Makes boxplots

exo_sets1_df = lapply(binsizes,FUN = logcount_data.frame,exo_sets1,mc.cores)
exo_sets1_df = mergeDataFrames(exo_sets1_df)
p1 = logcount_boxplot(exo_sets1_df)
pdf(file = file.path(figsdir,"chip_exo_logcounts_1300.pdf"),width = 9,height = 5)
print(p1);dev.off()

exo_sets2_df = lapply(binsizes,FUN = logcount_data.frame,exo_sets2,mc.cores)
exo_sets2_df = mergeDataFrames(exo_sets2_df)
p2 = logcount_boxplot(exo_sets2_df)
pdf(file = file.path(figsdir,"chip_exo_logcounts_900.pdf"),width = 9,height =5)
print(p2);dev.off()

pet_df = lapply(binsizes,FUN = logcount_data.frame,pet,mc.cores)
pet_df = mergeDataFrames(pet_df)
p3 = logcount_boxplot(pet_df)
pdf(file = file.path(figsdir,"chip_seq_pet_logcounts.pdf"),width = 9,height =5)
print(p3);dev.off()

set_df = lapply(binsizes,FUN = logcount_data.frame,set,mc.cores)
set_df = mergeDataFrames(set_df)
p4 = logcount_boxplot(set_df)
pdf(file = file.path(figsdir,"chip_seq_set_logcounts.pdf"),width = 4,height =5)
print(p4);dev.off()

# Generate sample summary info with depths updated by bam

sample.list = lapply(levels(sample.info$seq),function(x,samp)
  subset(samp,samp$seq==x),sample.info)
exo.depth = mclapply(exo,length,mc.cores=mc.cores)
pet.depth = mclapply(pet,length,mc.cores=mc.cores)
set.depth = mclapply(set,length,mc.cores=mc.cores)
depth = list(exo.depth,pet.depth,set.depth)
sample.list = lapply(1:3,function(i,sample.list,depth)assignDepth(sample.list[[i]],depth[[i]]),sample.list,depth)
names(sample.list) = c("Exo","PET","SET")

WriteXLS("sample.list",ExcelFileName = "inst/alignment/sample.summary_updated_with_bam.xls",
         SheetNames = names(sample.list))
