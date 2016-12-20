
rm(list = ls())

library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(parallel)
library(GenomicAlignments)
library(readr)

devtools:::load_all("~/Desktop/Docs/Code/mosaics")

figs_dir = "/p/keles/ChIPexo/volume3/ChIPexo/figs/NAR_review/MCF7_GC_mapp"

work_dir = "/p/keles/ChIPexo/volume4/carroll_data/human"

infiles = list.files(work_dir,full.names = TRUE,recursive = TRUE)
infiles = infiles[grep("sort",infiles)]
infiles = infiles[grep("bai",infiles,invert = TRUE)]

reads = mclapply(infiles,ChIPUtils:::create_reads,mc.cores =8)

chrom.sizes = read_delim("/p/keles/SOFTWARE/hg19.chrom.sizes",delim = "\t",col_names =FALSE)
chrom.sizes = GRanges(seqnames = chrom.sizes$X1,
                      ranges = IRanges(start = 1 , end = chrom.sizes$X2))

bins = mclapply(reads,function(x){
    ChIPUtils:::create_bins(200,reads = x,chrom = chrom.sizes, frag_len = 200)},mc.cores = 8)

binsDF = lapply(bins,function(x)tibble(chr = as.character(seqnames(x)),start = as.integer(start(x)),
                                       counts = as.integer(x$tagCounts)))

depth = lapply(reads,function(x)x@nReads)

write_bins <- function(bin,depth,file)
{
    line = paste(sep = " ","#","sequencing depth:",depth)
    write(line,file = file,append = TRUE)

    write_delim(bin,path = file,delim = "\t",append = TRUE)

}

u = mcmapply(write_bins,binsDF,depth,file.path(work_dir,"bins",
         paste0(basename(infiles),"_fragL",200,"_bin",200,".txt")),
         mc.cores = 8)

## first step, build data

read_bins <- function(file)
{
    Mfile  = file.path(work_dir,"extra","all_map_fragL200_bin200.txt")
    GCfile = file.path(work_dir,"extra","all_GC_fragL200_bin200.txt")
    Nfile = file.path(work_dir,"extra","all_N_fragL200_bin200.txt")
    readBins(c("chip","M","GC","N"),
             c(file,Mfile,GCfile,Nfile))
}

binfiles = list.files(file.path(work_dir,"bins"),full.names = TRUE)

bins = mclapply(binfiles,read_bins,mc.cores = 8)

## bin.exo <- readBins( c("chip","M","GC","N"),
## 	c( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/ChIP-exo/CTCF.bowtie_fragL200_bin200.txt",
## 	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/map50/map50_fragL200_bin200.txt",
## 	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/GC/GC_fragL200_bin200.txt",
## 	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/N/N_fragL200_bin200.txt" ) )

## bin.seq <- readBins( c("chip","input","M","GC","N"),
## 	c( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/ChIP-seq_Crawford/CTCF_Crawford.bowtie_fragL200_bin200.txt",
## 	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/Input_Crawford/Input_Crawford.bowtie_fragL200_bin200.txt",
## 	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/map50/map50_fragL200_bin200.txt",
## 	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/GC/GC_fragL200_bin200.txt",
## 	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/N/N_fragL200_bin200.txt" ) )


## Rlist <- system( "ls /u/w/e/welch/Desktop/Docs/packages/mosaics/R/*.R", intern=TRUE )
## for ( i in 1:length(Rlist) ) {
## 	source( Rlist[i] )
## }


calculate_stats <- function(bin)
{
    statM = mosaics:::.computeStat( Y=tagCount(bin), S=mappability(bin) )
    statGC = mosaics:::.computeStat( Y=tagCount(bin), S=gcContent(bin) )
    out = list(map = statM , GC = statGC)
    out
}

stats = mclapply(bins,calculate_stats,mc.cores =8 )

as_tibble <- function(stat)
{
    map = stat$map %>% as.data.frame %>% as.tbl %>%
        mutate(lab = "map") %>%
        arrange(desc(varYall))
    gc = stat$GC %>% as.data.frame %>% as.tbl %>%
        mutate(lab = "GC")
    rbind(map,gc) %>%
        mutate(lb = meanYall - 1.96 * sqrt(varYall / nitem),
               ub = meanYall + 1.96 * sqrt(varYall / nitem))

}


stats = lapply(stats,as_tibble)

theme_set(theme_bw())

generate_plot <- function(stat,name,ll)
{    
    p = ggplot(stat %>% filter(lab == ll & abs(ub - lb) < 10),
               aes(uS,ymin = lb,ymax = ub)) + geom_linerange()+
        geom_point(aes(x = uS,y = meanYall))+
        xlim(0,1)+
        ylab("Mean ChIP tag count")+ggtitle(name)
    if(ll == "GC"){
        p + xlab("GC content score")+
            geom_vline(xintercept = .5,colour = "red",linetype = 2)
    }else{
        p + xlab("Mappability score")+ylim(0,3)
    }
}

nms = c("ChIP-exo Rep-1","ChIP-exo Rep-2","ChIP-exo Rep-3",
        "ChIP-seq Rep-1","ChIP-seq Rep-2","ChIP-seq Rep-3",
        "ChIP-seq Input")



gc = mapply(generate_plot,stats,nms,MoreArgs = list("GC"),SIMPLIFY = FALSE)
map = mapply(generate_plot,stats,nms,MoreArgs = list("map"),SIMPLIFY = FALSE)


pdf("figs/NAR_review/Eukaryotic_biases/MCF7_GC.pdf")
u = lapply(gc,print)
dev.off()

pdf("figs/NAR_review/Eukaryotic_biases/MCF7_map.pdf")
u = lapply(map,print)
dev.off()

