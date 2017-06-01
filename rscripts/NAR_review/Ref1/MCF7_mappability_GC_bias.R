
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

infiles = infiles[grep("ERR",infiles)]
infiles = infiles[grep("txt",infiles,invert = TRUE)]

reads = mclapply(infiles,readGAlignments,param = NULL)
reads = lapply(reads,as,"GRanges")
reads = lapply(reads,resize, 200)

chrom.sizes = read_delim("/p/keles/SOFTWARE/hg19.chrom.sizes",
                         delim = "\t",col_names =FALSE)

bins = ChIPUtils::create_bins(200,chrom.sizes)

bins = lapply(reads,function(x){
    mcols(bins)$counts = countOverlaps(bins,x)
    bins})

binsDF = lapply(bins,
                function(x)
                    tibble(chr = as.character(seqnames(x)),
                           start = as.integer(start(x)),
                           counts = as.integer(x$counts)))

depth = lapply(reads,length)

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
binfiles = binfiles[grep("ERR",binfiles)]

bins = mclapply(binfiles,read_bins,mc.cores = 8)


computeStat <- function( Y, S )
{    
    # construct unique M/GC pair
    
    tmp <- cbind( Y, S )
    tempdf <- within( as.data.frame(tmp), {S <- factor(S)} )    
    ySubList <- with( tempdf, split(Y, S) )  
    
    # calculate strata-specific parameters
    
    yStatList <- lapply( ySubList, function(Y) {
            Y <- as.numeric( as.vector(Y) )
            nitem <- length(Y)
            nitemgeq0 <- length( Y[Y!=0] )
            p0 <- length(which( Y==0 )) / length(Y)
            p1 <- length(which( Y==1 )) / length(Y)
            meanYgeq0 <- mean( Y[Y!=0],na.rm = TRUE )
            varYgeq0 <- var( Y[Y!=0] ,na.rm = TRUE)
            medYgeq0 <- median( Y[Y!=0],na.rm = TRUE )
            meanYall <- mean(Y,na.rm = TRUE)
            varYall <- var(Y,na.rm = TRUE)
            
            result <- c( nitem, nitemgeq0, p0, p1,
                meanYgeq0, varYgeq0, medYgeq0, meanYall, varYall )
            
            return( result )
        }
    )
    yStatMat <- matrix( unlist(yStatList), ncol = 9, byrow = TRUE )
    
    # construct summary
    
    yStats <- list( uS = as.numeric(as.vector(names(ySubList))),
        nitem = as.numeric(as.vector(yStatMat[, 1])),
        nitemgeq0 = as.numeric(as.vector(yStatMat[, 2])),
        p0 = as.numeric(as.vector(yStatMat[, 3])),
        p1 = as.numeric(as.vector(yStatMat[, 4])), 
        meanYgeq0 = as.numeric(as.vector(yStatMat[, 5])), 
        varYgeq0 = as.numeric(as.vector(yStatMat[, 6])), 
        medYgeq0 = as.numeric(as.vector(yStatMat[, 7])),
        meanYall = as.numeric(as.vector(yStatMat[, 8])),
        varYall = as.numeric(as.vector(yStatMat[, 9]))
    )
        
    return(yStats)
}

calculate_stats <- function(bin)
{
    statM = computeStat( Y=tagCount(bin), S=mappability(bin) )
    statGC = computeStat( Y=tagCount(bin), S=gcContent(bin) )
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
        xlim(0,1)+theme_bw()+
        ylab("Mean ChIP tag count")+ggtitle(name)
    if(ll == "GC"){
        p + xlab("GC content score")+
            geom_vline(xintercept = .5,colour = "red",linetype = 2)
    }else{
        p + xlab("Mappability score")+ylim(0,3)
    }
}

nms = c("ChIP-exo Rep-1","ChIP-exo Rep-2","ChIP-exo Rep-3")

## ,
##         "ChIP-seq Rep-1","ChIP-seq Rep-2","ChIP-seq Rep-3",
##         "ChIP-seq Input")



gc = mapply(generate_plot,stats,nms,MoreArgs = list("GC"),SIMPLIFY = FALSE)
map = mapply(generate_plot,stats,nms,MoreArgs = list("map"),SIMPLIFY = FALSE)


pdf("figs/NAR_review/Eukaryotic_biases/MCF7_GC_exo.pdf")
u = lapply(gc,print)
dev.off()

pdf("figs/NAR_review/Eukaryotic_biases/MCF7_map_exo.pdf")
u = lapply(map,print)
dev.off()

