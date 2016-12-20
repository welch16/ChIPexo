

rm(list = ls())

library(ChIPexoQual)
library(parallel)
library(GenomicAlignments)

dr = "/p/keles/ChIPexo/volume4/exo_histone_data/BAM"
files = list.files(dr,full.names = TRUE)
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("histone1",files)]

cores = 20


##
reads = mclapply(files,readGAlignments,param = NULL,mc.cores = cores)
reads = mclapply(reads,as,"GRanges",mc.cores = cores)

exo = lapply(reads,function(x)ExoData(reads = x,mc.cores = cores))
names(exo) = paste0("Rep-",seq_along(exo))


library(ggplot2)
figsdr = "figs/NAR_review/exo_histone"

pdf(file.path(figsdr,"ARC_vURC.pdf"),width = 10,height = 5)
ARCvURCplot(exo)+xlim(0,5)
ARCvURCplot(exo,both.strand = TRUE)+xlim(0,5)
dev.off()

pdf(file.path(figsdr,"FSRdist.pdf"),height = 8)
FSRDistplot(exo,quantiles = c(.1,.25,.5,.75,.9),depth.values = seq_len(150))
FSRDistplot(exo,quantiles = c(.1,.25,.5,.75,.9),depth.values = seq_len(150),both.strand = TRUE)
dev.off()

pdf(file.path(figsdr,"RegionComp.pdf"),height = 8)
regionCompplot(exo,depth.values = seq_len(50))
dev.off()

pdf(file.path(figsdr,"beta_scores.pdf"))
paramDistBoxplot(exo,which.param = "beta1")
paramDistBoxplot(exo,which.param = "beta2")
dev.off()

## > sapply(exo,nreads)
##    Rep-1    Rep-2    Rep-3    Rep-4 
## 35,951,922 32,568,539 21,600,382 11,030,924 

library(ChIPUtils)

reads2 = mclapply(files,create_reads,mc.cores = cores)

pbc = mclapply(reads2,PBC,mc.cores = 4)
## > pbc %>% unlist
## [1] 0.4435131 0.3902492 0.4774308 0.5381422

genome = readDNAStringSet("/p/keles/ChIPexo/volume4/exo_histone_data/index/FASTA/yeast_S288C.fa")

library(data.table)
chr = seqnames(reads[[1]]) %>% as.character %>% unique

sizes = data.table(chr[1:16],
                   width(genome))


scc = lapply(reads2,strand_cross_corr,seq_len(250),sizes,TRUE)
scc = mapply(function(x,y)x[,rep := y],scc,paste0("Rep-",1:4),SIMPLIFY = FALSE) %>% rbindlist


pdf(file = file.path(figsdr,"SCC.pdf"),height = 5)
ggplot(scc,aes(shift,cross.corr,colour = rep))+geom_line()+
    scale_color_brewer(palette = "Set1")
dev.off()

## scc[,max(cross.corr) / min(cross.corr),by = rep]$V1
## > scc[,max(cross.corr) / min(cross.corr),by = rep]$V1
## [1]  5.650432  3.904779  9.174575 13.138201
