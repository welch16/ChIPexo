
rm(list = ls())

figsdr = "figs/NAR_review/palson"
dr = "/p/keles/ChIPexo/volume4/palsson_data/BAM"


library(ChIPexoQual)

files = list.files(dr,full.names = TRUE)
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("sort",files)]

library(parallel)

options(mc.cores = 20)

reads = mclapply(files,readGAlignments,param = NULL)

library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

reads = reads %>% mclapply(as,"GRanges")

exo = lapply(reads,function(x)ExoData(reads = x,nregions = 200))

names(exo) = files %>% sapply(function(x)gsub(".sort.bam","",basename(x))) %>%
    sapply(function(x)gsub("exo_K12_MG1655_","",x))


DT = tibble(ff = names(exo),
         nreads  = reads %>% sapply(length),
         nregions = exo %>% sapply(length) ,
         beta1_median = exo %>% sapply(function(x)median(paramDist(x)$beta1)),
         beta1_mean = exo %>% sapply(function(x)mean(paramDist(x)$beta1)),
         beta2_median = exo %>% sapply(function(x)-median(paramDist(x)$beta2)),
         beta2_mean = exo %>% sapply(function(x)-mean(paramDist(x)$beta2)))

pdf(file.path(figsdr,"Ecoli_Nreads_vs_nregions.pdf"))
print( DT
  %>%
    ggplot(aes(nreads,nregions))+geom_point()+geom_text_repel(aes(label = ff)) +
    scale_y_log10()
)
dev.off()



pdf(file.path(figsdr,"Ecoli_palson_ARC_vURC.pdf"),height = 10,width = 12)
ARCvURCplot(exo)
dev.off()


pdf(file.path(figsdr,"Ecoli_palson_FSRdist.pdf"),height = 14)
FSRDistplot(exo,quantiles = c(.1,.25,.5,.75,.9),depth.values = seq_len(50))
dev.off()


pdf(file.path(figsdr,"Ecoli_palson_RegionComp.pdf"),height = 14)
regionCompplot(exo,depth.values = seq_len(50))
dev.off()

pdf(file.path(figsdr,"Ecoli_palson_qcscores.pdf"),width = 10)
paramDistBoxplot(exo)
dev.off()

## # A tibble: 14 × 7
##           ff  nreads nregions beta1_median beta1_mean beta2_median beta2_mean
##        <chr>   <int>    <int>        <dbl>      <dbl>        <dbl>      <dbl>

## 1  gadE_rep1 2 473 587      221     2.157916   2.154317   0.21455042 0.21341076
## 2  gadE_rep2 2 256 541      216     1.995088   1.991756   0.16153837 0.16039322
## 3  gadW_rep1 2 848 471      368     4.562672   4.477912   0.24763411 0.23157096
## 4  gadW_rep2 2 490 801      247     3.097940   3.073969   0.20267380 0.19718388
## 5  gadX_rep1 2 268 259      397    10.587929   8.419101   1.67951707 1.22908339
## 6  gadX_rep2 2 489 838     1043     6.378889   9.354647   0.41822388 0.88267107

## 7  oxyr_rep1 2 174 518     1145     2.415488   2.380502   0.18416876 0.17721439
## 8  oxyr_rep2 3 319 297     1671     3.820536   3.966763   0.20999988 0.24970408
## 9  rpoS_rep1 5 396 499     4313     7.392655   8.515317   1.06991209 1.37657525
## 10 rpoS_rep2 5 154 827     5712     6.863224   7.965851   0.98887748 1.38058384
## 11 soxr_rep1 1 494 706    43974     5.604541   5.642230   0.03352554 0.04474422
## 12 soxr_rep2 1 567 538    45893     6.628990   7.634576   0.02805468 0.13164397
## 13 soxs_rep1 2 088 467    38777     6.577314   6.937644   0.05803393 0.10233846
## 14 soxs_rep2 2 283 272    35406     6.307526   6.939170   0.07059690 0.15787385


## it seems that the GADEWX samples have a low ammount of peaks, hence it is reasonable to assume that
## the reads are located around those regions and a bit more

library(ChIPUtils)

reads_cu = files %>% mclapply(create_reads)


pbc = reads_cu %>% mclapply(PBC) %>% unlist

sizes = read.table("/p/keles/ChIPexo/volume7/Landick/K12/K12_size") %>% data.table


scc <- lapply(reads_cu, strand_cross_corr,shift = seq_len(300),
   chrom.sizes = sizes,parallel = FALSE)
 

scc = mapply(function(x,y)x[,sample := y],scc,DT$ff,SIMPLIFY = FALSE)
scc = rbindlist(scc)

scc = scc %>% as.tbl

library(ggplot2)
library(tidyr)


pdf(file.path(figsdr,"Ecoli_SCC.pdf"),height = 12,width = 6)
scc %>% separate(sample,into = c("TF","Rep"),sep = "_") %>%    
    ggplot(aes(shift,cross.corr,colour = Rep))+geom_line()+
      facet_grid(TF ~ . , scales = "free_y")+scale_color_brewer(palette = "Dark2")
dev.off()

