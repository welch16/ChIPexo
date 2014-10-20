

library(GenomicAlignments)
library(parallel)
filesdir = "/p/keles/ChIPexo/volume3/RenData/BAMfiles"


files = list.files(filesdir)
files = files[grep("bam",files)]
files = files[grep("sort",files)]
files = files[!grepl("bai",files)]

files1 = files[grep("AY552",files)]
files2 = files[grep("AY553",files)]
files3 = files[grep("AY554",files)]


reads = mclapply(file.path(filesdir,files2),FUN = readGAlignmentsFromBam,param = NULL,mc.cores = 2)
reads2 = mclapply(reads,function(x)as(x,"GRanges"),mc.cores=2)
names(reads2) = files2

reads = mclapply(file.path(filesdir,files1),FUN = readGAlignmentsFromBam,param = NULL,mc.cores = 2)
reads1 = mclapply(reads,function(x)as(x,"GRanges"),mc.cores=2)
names(reads1) = files1

reads = mclapply(file.path(filesdir,files3),FUN = readGAlignmentsFromBam,param = NULL,mc.cores = 2)
reads3 = mclapply(reads,function(x)as(x,"GRanges"),mc.cores=2)
names(reads3) = files3


datadir = "data"
save(file = file.path(datadir,"Ren_reads.RData"),list = c("reads1","reads2","reads3"))

