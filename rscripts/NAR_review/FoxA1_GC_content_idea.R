
rm(list = ls())

library(ChIPexoQual)
library(GenomicAlignments)
library(BSgenome.Mmusculus.UCSC.mm9)

in_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files <- list.files(in_dir,full.names = TRUE)
files <- files[grepl("sort",files) & !grepl("bai",files)]

library(parallel)

options(mc.cores = 20)

reads <- mclapply(files,readGAlignments,param = ScanBamParam(what = "seq"))

sequences <- lapply(reads,function(x)mcols(x)$seq)

library(GenomicRanges)

cleanReads <- function(read)
{
    read = as(read,"GRanges")
    mcols(read) = NULL
    read
}

clean_reads <- mclapply(reads,cleanReads)

library(BiocParallel)

exo_regions <- lapply(clean_reads,function(y)ExoData(reads = y,mc.cores = 22))

overlaps <- mcmapply(findOverlaps,clean_reads,exo_regions,SIMPLIFY = FALSE)

library(Biostrings)

freqs <- mclapply(sequences,oligonucleotideFrequency,width = 1)

gc_content_from_seq <- function(freq)
{
  (freq[,"C"] + freq[,"G"]) / rowSums(freq)
}

gc_content <- mclapply(freqs,gc_content_from_seq)

library(dplyr)
library(magrittr)

overlap_gc_tibble <- function(overlap,gc)
{
    tibble( queryHits = queryHits(overlap),subjectHits = subjectHits(overlap),gc)
}

overlaps <- mapply(overlap_gc_tibble,overlaps,gc_content,SIMPLIFY = FALSE)

gc <- mclapply(overlaps,
             function(x){
                 x %>% group_by(subjectHits) %>% summarize(gc = mean(gc))
             })

add_gc <- function(exo_reg,gc0)
{

    mcols(exo_reg)$gc = gc0$gc
    exo_reg
}

exo_regions <- mapply(add_gc,exo_regions,gc,SIMPLIFY = FALSE)
names(exo_regions) <- c("Rep3","Rep1","Rep2")

exo = mapply(function(a,b){
    a = mcols(a) %>% as.data.frame %>% as.tbl %>%
        mutate(rep = b)},exo_regions,names(exo_regions),SIMPLIFY = FALSE)

exo = bind_rows(exo)

library(scales)
library(ggplot2)
library(viridis)
library(hexbin)



scatter_plot <- function(exo,x,y)
{
    r = viridis(1e3 , option = "D")
  ggplot(exo,aes_string(x,y))+stat_binhex(bins = 70)+
                scale_fill_gradientn(colours = r,trans = 'log10',
                                     labels=trans_format('log10',math_format(10^.x)) )+
                            theme_bw()+theme(legend.position = "top")+
      facet_grid( . ~ rep)+
      geom_abline(slope = 0,intercept = .5,colour = "orange",linetype = 2)
}

figs_dir <- "figs/NAR_review/"

pdf(file.path(figs_dir,"FoxA1_GC_test.pdf"),width = 9,height = 4)
scatter_plot(exo,"ARC","gc")+xlim(0,4)
scatter_plot(exo,"URC","gc")
scatter_plot(exo,"FSR","gc")
dev.off()
