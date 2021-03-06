rm(list = ls())

library(ChIPexoQual)
library(magrittr)

dr = "/p/keles/ChIPexo/volume4"

files = list.files(dr,recursive = TRUE,full.names = TRUE)
files = files[grep("carroll",files)]
files = files[grep("bam",files)]
files = files[grep("mouse",files)]
files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("txt",files,invert = TRUE)]

library(parallel)

options("mc.cores" = 22)

reads = files %>% mclapply(readGAlignments,param = NULL)
reads = reads %>% mclapply(as,"GRanges")
names(reads) = paste0("Rep",c(3,1,2))

exo = lapply(reads,function(x)ExoData(reads =  x))

peakfiles = list.files(dr,full.names = TRUE,recursive = TRUE,pattern = "peaks")
peakfiles = peakfiles[grep("mouse",peakfiles)]

library(readr)
library(dplyr)
library(data.table)

peaks = mclapply(peakfiles,read_delim,delim = " ",col_names = FALSE)
peaks = lapply(peaks,function(x){
    x = x %>% select(X1,X2,X3)
    setnames(x,names(x),c("seqnames","start","end"))
    x = as.data.table(x)
  return(ChIPUtils::dt2gr(x))})
names(peaks) = paste0("Rep",c(3,1,2))

## peak columns:
## chrID peakStart peakStop peakSize logAveP logMinP aveLogP aveChipCount maxChipCount map GC 

## join all peak regions together
all_peaks = Reduce(c,peaks) %>% reduce

exo_peaks = mclapply(exo,subsetByOverlaps,all_peaks)
readlength = reads %>% sapply(function(x)x %>% width %>% median)

## common sense filter

## remove chrM
exo_peaks = exo_peaks %>% mclapply(function(x)x[as.character(seqnames(x)) != "chrM"])

## width analysis

exo_peaks %>% lapply(function(x)width(x) %>% summary)

## $Rep3
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   34.00   36.00   53.00   80.97   98.00 3399.00 

## $Rep1
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    35.0    60.0   233.0   254.8   400.0  3390.0 

## $Rep2
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    35.0    39.0    86.0   147.8   234.0  3388.0 

rl = median(readlength)

K = 3
exo_peaks = exo_peaks %>% lapply(function(x)x[width(x) >= K * rl])

## load fimo files
fimofiles = list.files(dr,pattern = "fimo",recursive = TRUE,full.names = TRUE)
fimofiles = fimofiles[grep("carroll",fimofiles)]
fimofiles = fimofiles[grep("txt",fimofiles)]
fimofiles = fimofiles[grep("FoxA1",fimofiles)]
fimo = lapply(fimofiles,read_delim,delim = "\t")
names(fimo) = paste0("Rep",seq_along(fimo))

library(ggplot2)
library(tidyr)

## scores boxplots
fimo = mapply(function(x,y)x %>% mutate(replicate = y),fimo,names(fimo),SIMPLIFY = FALSE)
all_fimo = fimo %>% bind_rows

theme_set(theme_bw())

score_boxplot <- function(all_fimo,topKcuts,ylims = NULL)
{
  
    scores = lapply(topKcuts,function(x){
        all_fimo %>% arrange_( ~ desc(score)) %>% group_by_(~replicate ) %>%
            top_n(n = x,wt = score) %>% mutate( K = x)}) %>% bind_rows
    plot = scores %>%
        ggplot(aes(replicate , score, fill = replicate))+ geom_boxplot()+
        scale_fill_brewer(palette = "Pastel1",name = "Sample") +
        theme(legend.position = "top",
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank())+
        facet_grid( . ~ K )
    if(!is.null(ylims))plot = plot + ylim(ylims)
    plot
        
        

}

pdf("figs/NAR_review/FoxA1_fimo_analysis_around_peaks_full.pdf",width = 9,height = 6)
score_boxplot(all_fimo,c(100,500,1e3,2e3), ylims = c(12.5,16.5)) ; dev.off()

## profiles

maxQval = 0.05
window_length = 20
sm = 1

## reverse_cover = reads %>% mclapply(function(x){
##     subset(x,as.character(strand(x)) == "-") %>% coverage } )

## forward_cover = reads %>% mclapply(function(x){
##     subset(x,as.character(strand(x)) == "+") %>% coverage } )

profile_table <- function(fimo_repl,reads,repl,maxQval,windowLength,smooth)
{
    bwd = subset(reads,as.character(strand) == "-")
    bwd = resize(bwd,width = 1) %>% coverage
    fwd = subset(reads,as.character(strand) == "+")
    fwd = resize(fwd,width = 1) %>% coverage
    depth = reads %>% length
    
    regions = fimo_repl %>% filter(`q-value` <= maxQval) %>%
        separate(`sequence name`,into = c("chr","rStart"),sep = "\\:") %>%
        separate(rStart , into = c("rstart","rend"),sep = "-") %>%
        mutate(rstart =as.numeric(rstart) , rend = as.numeric(rend)) 

    fwdranges = regions %>% filter(strand == "+") %>%
        mutate(start = start + rstart, stop  = start ) %>% DataFrame %>%
        makeGRangesFromDataFrame
    revranges = regions %>% filter(strand == "-") %>%
        mutate(start = stop + rstart, stop  = start ) %>% DataFrame %>%
        makeGRangesFromDataFrame

    fwdranges = resize(fwdranges,2*window_length + 1,fix = "center")
    revranges = resize(revranges,2*window_length + 1,fix = "center")
    
    coord = seq( - window_length, window_length)
    fwdcounts = fwd[fwdranges] %>% lapply(as.vector)
    bwdcounts = bwd[revranges] %>% lapply(as.vector)

    fwdcounts = fwdcounts %>% lapply(function(x)tibble(coord,counts = x)) %>%
        bind_rows %>% mutate(strand = "+")
    bwdcounts = bwdcounts %>% lapply(function(x)tibble(coord,counts = x)) %>%
        bind_rows %>% mutate(strand = "-")

    rbind(fwdcounts,bwdcounts) %>% mutate(replicate = repl,signal = 1e9 * counts / depth)
}





profiles = mcmapply(profile_table,fimo,reads,names(reads),
       maxQval,window_length,sm ,SIMPLIFY = FALSE) %>% bind_rows


## profile_table_strand  <- function(fimo_repl,reads,repl,maxQval,windowLength,smooth,strand)
## {
##     browser()
##     bwd = subset(reads,as.character(strand) == "-")
##     bwd = resize(bwd,width = 1) %>% coverage
##     fwd = subset(reads,as.character(strand) == "+")
##     fwd = resize(fwd,width = 1) %>% coverage
##     depth = reads %>% length
    
##     regions = fimo_repl %>% filter(`q-value` <= maxQval) %>%
##         separate(`sequence name`,into = c("chr","rStart"),sep = "\\:") %>%
##         separate(rStart , into = c("rstart","rend"),sep = "-") %>%
##         mutate(rstart =as.numeric(rstart) , rend = as.numeric(rend)) 

##     fwdranges = regions %>% filter(strand == "+") %>%
##         mutate(start = start + rstart, stop  = start ) %>% DataFrame %>%
##         makeGRangesFromDataFrame
##     revranges = regions %>% filter(strand == "-") %>%
##         mutate(start = stop + rstart, stop  = start ) %>% DataFrame %>%
##         makeGRangesFromDataFrame

##     fwdranges = resize(fwdranges,2*window_length + 1,fix = "center")
##     revranges = resize(revranges,2*window_length + 1,fix = "center")
    
##     coord = seq( - window_length, window_length)
##     fwdcounts = fwd[fwdranges] %>% lapply(as.vector)
##     bwdcounts = bwd[revranges] %>% lapply(as.vector)

##     fwdcounts = fwdcounts %>% lapply(function(x)tibble(coord,counts = x)) %>%
##         bind_rows %>% mutate(strand = "+")
##     bwdcounts = bwdcounts %>% lapply(function(x)tibble(coord,counts = x)) %>%
##         bind_rows %>% mutate(strand = "-")

##     rbind(fwdcounts,bwdcounts) %>% mutate(replicate = repl,signal = 1e9 * counts / depth)
## }


## profiles_strand  = mapply(profile_table_strand,fimo,reads,names(reads),
##        maxQval,window_length,sm,"+" ,SIMPLIFY = FALSE) %>% bind_rows


strand_imbalance_fimo <- function(exo,fimo,repl)
{

    exo_tbl = exo %>% as.data.frame %>% as.tbl %>%
        mutate(match = paste0(seqnames,":",start,"-",end)) %>%
        select(match , everything())
    fimo = fimo %>% rename(match = `sequence name`)

    fimo_summary = fimo %>% group_by(match) %>%
        summarize(
            nmotif = n(),
            minPval = min(`p-value`),
            maxScore = max(score)
        )

    exo_tbl %>% left_join(fimo_summary,by ="match") %>%
        mutate(nmotif = ifelse(is.na(nmotif),0,nmotif),repl)

  
}


imbalance = mapply(strand_imbalance_fimo,exo_peaks,fimo,names(fimo),SIMPLIFY = FALSE) %>% bind_rows


r4 = c("darkblue","firebrick3","blue","red","lightblue","lightpink")

pal = RColorBrewer::brewer.pal(n = 8,"Set2")

pdf("figs/NAR_review/FoxA1_fimo_analysis_around_peaks.pdf",width = 9,height = 6)
print(score_boxplot(all_fimo,c(50,100,250,500,1e3,2e3), ylims = c(12.5,16.5)))
print(profiles %>% group_by(coord,replicate,strand) %>% summarize(signal = mean(counts)) %>%
    mutate(class = paste(replicate,strand,sep = ":")) %>%
    ggplot(aes(coord,signal,colour = class))+
    geom_line(size = 1.1)+
    scale_color_manual(values = r4)+
    theme(legend.position = "top")+
  xlab("Position around motif start")+ylab("Average counts"))
print(imbalance %>% mutate(nmotif = forcats::fct_collapse(as.factor(nmotif),">=3" = as.character(seq(3,100)))) %>%
    ggplot(aes(FSR,fill = nmotif))+
    geom_histogram(aes(y = ..count..),colour = "black",bins = 25) + facet_grid( . ~ repl)+
    theme(legend.position = "top")+
    scale_fill_brewer(palette  = "Set2",name = "Number of motifs matches in region")+
    ylab("ChIP-exo regions overlapping with peaks")+xlab("Forward Strand Ratio")+
    geom_vline(xintercept = .5,linetype = 2,colour = "darkgrey"))
print(imbalance %>% filter(nmotif > 0) %>%
    mutate(nmotif = forcats::fct_collapse(as.factor(nmotif),">=3" = as.character(seq(3,100)))) %>%
    ggplot(aes(FSR,fill = nmotif))+
    geom_histogram(aes(y = ..count..),colour = "black",bins = 25) + facet_grid( . ~ repl)+
    theme(legend.position = "top")+
    scale_fill_manual(values = pal[-1],name = "Number of motifs matches in region")+
    ylab("ChIP-exo regions overlapping with peaks")+xlab("Forward Strand Ratio")+
    geom_vline(xintercept = .5,linetype = 2,colour = "darkgrey"))
print(imbalance %>% filter(nmotif > 1) %>%
    mutate(nmotif = forcats::fct_collapse(as.factor(nmotif),">=3" = as.character(seq(3,100)))) %>%
    ggplot(aes(FSR,fill = nmotif))+
    geom_histogram(aes(y = ..count..),colour = "black",bins = 25) + facet_grid( . ~ repl)+
    theme(legend.position = "top")+
    scale_fill_manual(values = pal[-(1:2)],name = "Number of motifs matches in region")+
    ylab("ChIP-exo regions overlapping with peaks")+xlab("Forward Strand Ratio")+
    geom_vline(xintercept = .5,linetype = 2,colour = "darkgrey"))
dev.off()

