
rm(list = ls())

library(ChIPexoQual)
library(parallel)
library(GenomicAlignments)

dr = "/p/keles/ChIPexo/volume4/exo_histone_data/BAM"
files = list.files(dr,full.names = TRUE)
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("histone1",files)]

cores = 20

library(ChIPUtils)
library(magrittr)

reads = files %>% lapply(create_reads)

## Depth

depth = reads %>% lapply(nreads)
## > depth
## [[1]]
## [1] 35,951,922

## [[2]]
## [1] 32,568,539

## [[3]]
# [1] 21,600,382

## [[4]]
## [1] 11,030,924

## PBC
pbc = reads %>% sapply(PBC)

## > pbc
## [1] 0.4435131 0.3902492 0.4774308 0.5381422



## SSD
library(htSeqTools)

rr = reads %>% lapply(function(x)rbindlist(c(readsF(x),readsR(x)))) %>%
    lapply(function(x)GRanges(seqnames = x$seqnames, ranges = IRanges(start = x$start,end = x$end)))

ssd2 = rr %>% sapply(ssdCoverage)
gini2 = rr %>% sapply(giniCoverage,mc.cores = 22)

## > ssd2
## [1] 24.17524 17.27585 29.04963 30.90528
## > gini2
##                 [,1]      [,2]      [,3]      [,4]
## gini        0.619365 0.3776861 0.6194098 0.5477689
## gini.adjust 0.562344 0.3177126 0.5458642 0.4450202



## NSC and RSC

library(data.table)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

seqlengths(Scerevisiae) -> sl

sizes = data.table(V1 = rr[[1]] %>% seqlengths %>% names,V2 = sl)

scc = lapply(reads,strand_cross_corr,seq_len(300),sizes,TRUE)
names(scc) = paste0("Rep-",seq_along(scc))

scc = mapply(function(x,y)x[,rep := y],scc,names(scc),SIMPLIFY = FALSE) %>% rbindlist

figs = "figs/NAR_review"

pdf(file = file.path(figs,"exohistone_SCC.pdf"),height = 5)
ggplot(scc,aes(shift,cross.corr,colour = rep))+geom_line()+
    scale_color_brewer(palette = "Set1")
dev.off()



rl = rr %>% lapply(width) %>% sapply(median)




RSC1 <- function(scc,read_length){
  out <- scc[,max(cross.corr)] / scc[shift == read_length, (cross.corr)]
  return(out)
}

RSC2 <- function(scc,read_length){
  mm <- scc[,min(cross.corr)]
  out <- (scc[,max(cross.corr)] - mm) /( scc[shift == read_length, (cross.corr)] - mm)
  return(out)
}

RSC1(scc[rep == "Rep-1"],rl[1])
RSC1(scc[rep == "Rep-2"],rl[2])
RSC1(scc[rep == "Rep-3"],rl[3])
RSC1(scc[rep == "Rep-4"],rl[4])

## > RSC1(scc[rep == "Rep-1"],rl[1])
## [1] 1.063164
## > RSC1(scc[rep == "Rep-2"],rl[2])
## [1] 1.211272
## > RSC1(scc[rep == "Rep-3"],rl[3])
## [1] 1.165552
## > RSC1(scc[rep == "Rep-4"],rl[4])
## [1] 1.113772


RSC2(scc[rep == "Rep-1"],rl[1])
RSC2(scc[rep == "Rep-2"],rl[2])
RSC2(scc[rep == "Rep-3"],rl[3])
RSC2(scc[rep == "Rep-3"],rl[4])

## RSC2(scc[rep == "Rep-1"],rl[1])
## [1] 1.074857
## > RSC2(scc[rep == "Rep-2"],rl[2])
## [1] 1.292136
## > RSC2(scc[rep == "Rep-3"],rl[3])
## [1] 1.183498
## > RSC2(scc[rep == "Rep-3"],rl[4])
## [1] 1.183498


scc[,max(cross.corr)/min(cross.corr),by = rep]
##      rep        V1
## 1: Rep-1  6.806251
## 2: Rep-2  4.375933
## 3: Rep-3 11.917300
## 4: Rep-4 14.923625


