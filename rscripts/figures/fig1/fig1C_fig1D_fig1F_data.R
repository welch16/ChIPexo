
rm( list = ls())

library(readr)
library(GenomicAlignments)
library(dplyr)
library(GenomicRanges)
library(parallel)

options(mc.cores = 10)

Rlist <- system( "ls /u/w/e/welch/Desktop/Docs/packages/mosaics/R/*.R", intern=TRUE )
for ( i in 1:length(Rlist) ) {
	source( Rlist[i] )
}

fragL = 200
binSize = 200

dr = "/p/keles/ChIPexo/volume4/pugh_data/mace"
files = list.files(dr,full.names = TRUE,pattern = "bam")
files = files[grep("bai",files,invert = TRUE)]
names(files) = paste0("Rep",seq_along(files))

reads = files %>% lapply(readGAlignments,param = NULL)
reads = reads %>% lapply(as,"GRanges")

reads = lapply(reads,function(x)resize(x,fragL))


chrom.sizes = "/p/keles/SOFTWARE/hg19.chrom.sizes"
Mfile = "/p/keles/ChIPexo/volume4/carroll_data/human/extra/hg19_mappability_fragL200_bin200_36mer_single/all_map_fragL200_bin200.txt"
GCfile = "/p/keles/ChIPexo/volume4/carroll_data/human/extra/hg19_GC_fragL200_bin200/all_GC_fragL200_bin200.txt"
Nfile = "/p/keles/ChIPexo/volume4/carroll_data/human/extra/hg19_N_fragL200_bin200/all_N_fragL200_bin200.txt"

sizes = read_delim(chrom.sizes,delim = "\t",col_names = FALSE)

read_file <- function(file)
{
    read_delim(file,delim = "\t",col_names = c("seqnames","start","val"),
               col_types = cols(seqnames = col_character(),
                                start = col_integer(),
                                val = col_double()))
}

mappp = read_file(Mfile)
gccontent = read_file(GCfile)
nDNA = read_file(Nfile)


mapbins = GRanges(mappp$seqnames,ranges = IRanges(start = mappp$start + 1, width = binSize))
mapbins$M = mappp$val

gcbins = GRanges(gccontent$seqnames,ranges = IRanges(start = gccontent$start + 1,width = binSize))
gcbins$GC = gccontent$val

nbins = GRanges(nDNA$seqnames,ranges = IRanges(start = nDNA$start + 1,width = binSize))
nbins$N = nDNA$val

sizes2 = GRanges(sizes$X1,IRanges(start = 1,width = sizes$X2))
bins = ChIPUtils::create_bins(binSize,chrom = sizes2)

bins$M = 0
bins$N = 0
bins$GC = 0

bins = shift(bins,1)

add_col <- function(which,bins,col)
{
    nn = nearest(bins,col)

    mcols(bins)[[which]] = mcols(col[nn])[[which]]

    bins
    
}

bins = add_col("M",bins,mapbins)
bins = add_col("GC",bins,gcbins)
bins = add_col("N",bins,nbins)

tagcounts = lapply(reads,function(x)countOverlaps(bins,x))

library(mosaics)

mbins = lapply(tagcounts,function(x,bins){
    new("BinData",
        chrID = as.character(seqnames(bins)),
        coord = end(bins),
        tagCount = x,
        mappability = bins$M,
        gcContent = bins$GC,
        input = bins$N,
        dataType = c("chip","M","GC","N"))},bins)




statM = mbins %>% lapply(function(x).computeStat( Y=tagCount(x), S=mappability(x) ))
statGC = mbins %>% lapply(function(x).computeStat( Y=tagCount(x), S=gcContent(x)))

M = statM %>% lapply(function(x){
    tibble(uS = x$uS,meanYall = x$meanYall,
           varYall = x$varYall, nitem = x$nitem)})
GC = statGC %>% lapply(function(x){
    tibble(uS = x$uS,meanYall = x$meanYall,
           varYall = x$varYall, nitem = x$nitem)})

M = mapply(function(x,y)x %>% mutate(sample = y),M,names(files),SIMPLIFY = FALSE) %>% bind_rows
GC = mapply(function(x,y)x %>% mutate(sample = y),GC,names(files),SIMPLIFY = FALSE) %>% bind_rows


quant = 1.96
M = M %>% select(sample,everything()) %>%
    mutate(lb = meanYall - quant * sqrt(varYall / nitem),
           ub = meanYall + quant * sqrt(varYall / nitem))

GC = GC %>% select(sample,everything()) %>%
    mutate(lb = meanYall - quant * sqrt(varYall / nitem),
           ub = meanYall + quant * sqrt(varYall / nitem))

write_tsv(M,"data/figures/fig1/fig1C_mappability.tsv")
write_tsv(GC,"data/figures/fig1/fig1D_GCcontent.tsv")




  ## lb = statM$meanYall-1.96*sqrt(statM$varYall/statM$nitem),
  ## ub = statM$meanYall+1.96*sqrt(statM$varYall/statM$nitem))

## pdf(file = "figs/for_paper/eukaryotic_bias_CTCF.pdf",width = 4 , height = 4)
## ggplot( M , aes( x = mapp ,ymin = lb , ymax = ub))+geom_linerange()+xlim(0,1)+
##   geom_point(aes( x = mapp, y = mean),size = 2)+theme_bw()+
##   theme(plot.title = element_text(hjust = 0))+ggtitle("C")+
##   coord_cartesian(xlim = c(0,1), ylim = c(-.1,7))+
##   xlab("Mappability score")+ylab("Mean ChIP tag count")
## ggplot( GC , aes( x = gc ,ymin = lb , ymax = ub))+geom_linerange()+xlim(0,1)+
##   geom_point(aes( x =gc, y = mean),size = 2)+theme_bw()+
##   theme(plot.title = element_text(hjust = 0))+ggtitle("D")+
##   coord_cartesian(xlim = c(0,1), ylim = c(-.2,15))+
##   xlab("GC content score")+ylab("Mean ChIP tag count")
## dev.off()




## pdf( paste(dir_out,"seq_bias_CTCF_mappability.pdf",sep="") )        
## plot( statM$uS, statM$meanYall,
##     xlab='Mappability score', ylab='Mean ChIP tag count', 
##     #main='Mappability score vs. Mean ChIP tag count',
##     ylim=quantile( statM$meanYall, prob=c(0.05,0.95) ),
##     cex.axis=1.5, cex.lab=1.5 )
## segments( statM$uS, statM$meanYall, 
##     statM$uS, statM$meanYall+1.96*sqrt(statM$varYall/statM$nitem) )
## segments( statM$uS, statM$meanYall, 
##     statM$uS, statM$meanYall-1.96*sqrt(statM$varYall/statM$nitem) )
## dev.off()

## pdf( paste(dir_out,"seq_bias_CTCF_GC.pdf",sep="") )
## plot( statGC$uS, statGC$meanYall,
##     xlab='GC content score', ylab='Mean ChIP tag count', 
##     #main='GC content score vs. Mean ChIP tag count',
##     ylim=quantile( statGC$meanYall, prob=c(0.05,0.95) ),
##     cex.axis=1.5, cex.lab=1.5 )
## segments( statGC$uS, statGC$meanYall, 
##     statGC$uS, statGC$meanYall+1.96*sqrt(statGC$varYall/statGC$nitem) )
## segments( statGC$uS, statGC$meanYall, 
##     statGC$uS, statGC$meanYall-1.96*sqrt(statGC$varYall/statGC$nitem) )
## dev.off()
## #plot( bin.seq, plotType="M" )
## #plot( bin.seq, plotType="GC" )
	
