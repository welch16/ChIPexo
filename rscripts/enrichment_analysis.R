 
rm(list = ls())

library(data.table)
library(GenomicRanges)
library(parallel)
library(ggplot2)

dir <- "/p/keles/ChIPexo/volume3/Analysis/Carroll/mouse/FoxA1-rep3/data"
figs_dir <- "figs/enrichment/Carroll/mouse/FoxA1-rep3"

build_files <- function(dir)
{
  files <- list.files(dir)
  patterns <- c("depth","reads_by_region","regions","summary_stat")
  files <- files[sapply(patterns,function(x,files)
                        grep(x,files),files)]
  return(files)
}

files <- build_files(dir)

mc <- detectCores()
## depth , reads , regions and summary_stats

load(file.path(dir,files[1]))  ## depth
load(file.path(dir,files[2]))  ## reads_table
load(file.path(dir,files[3]))  ## regions
load(file.path(dir,files[4]))  ## summary_stats

reads_table <- lapply(reads_table,
  function(x){
    x[,match := paste(seqnames,region,sep ="_")]
    x[,region := NULL]
    return(x)})

reads_table <- do.call(rbind,reads_table)
setkey(reads_table, match)

chr <- names(regions)
regions <- mclapply(chr,function(x,regions){
  out <- data.table(seqnames = x ,
                    start = start(regions[[x]]),
                    end = end(regions[[x]]))
  out[,match := paste(seqnames , 1:nrow(out),sep = "_")]
  return(out)},regions,mc.cores = mc)
regions <- do.call(rbind,regions)
setkey(regions,match)

summary_stats <- mclapply(summary_stats,
  function(x){
    x[,region := paste(chrID,region,sep = "_")]
    x[,chrID := NULL]
    nms <- names(x)
    setnames(x,nms,c("match",nms[-1]))
    return(x)},mc.cores = mc)
summary_stats <- do.call(rbind,summary_stats)
setkey(summary_stats,match)

stats <- merge(regions,summary_stats,by  = "match")

rm(regions,summary_stats)


dt2ir <- function(x)IRanges(start = x[,(start)],end = x[,(end)])


dt2gr <- function(x)
  GRanges(seqnames = x[,(seqnames)],
          ranges = dt2ir(x),
          strand = "*")

normalize.tagcounts <- function(counts,depth)return(counts*1e6 / depth)

Rle2dt <- function(rle_data)
{
  coord = cumsum(runLength(rle_data))
  counts = runValue(rle_data)
  return(data.table(coord=coord,counts = counts))
}
  
get_profile <- function(key,stats,reads_table,fl = 1,depth = 1)
{

  chr <- stats[key,(seqnames)]
  start <- stats[key,(start)]
  end <- stats[key,(end)]
  fwd <- coverage(resize(dt2ir(reads_table[key][strand == "+"]),fl))
  bwd <- coverage(resize(dt2ir(reads_table[key][strand == "-"]),fl))

  fwd <- Rle2dt(fwd)
  bwd <- Rle2dt(bwd)
  fwd[,strand := "+"]
  bwd[,strand := "-"]

  fwd[,counts := normalize.tagcounts(counts,depth)]
  bwd[,counts := -normalize.tagcounts(counts,depth)]

  M <- max(fwd[,max(counts)],-min(bwd[,counts]))

  fwd2 <- data.table(coord = start:end,counts = 0 ,strand = "+")
  bwd2 <- data.table(coord = start:end,counts = 0 ,strand = "-")

  if(nrow(fwd) > 0){
    fwd2[ between(coord,fwd[,min(coord)],fwd[,max(coord)]),
         counts := approxfun(x = fwd[,(coord)],y = fwd[,(counts)],method = "constant")(coord)]
  }
  if(nrow(bwd) > 0){
    bwd2[ between(coord,bwd[,min(coord)],bwd[,max(coord)]),
         counts := approxfun(x = bwd[,(coord)],y = bwd[,(counts)],method = "constant")(coord)]
  }   


  
  p = ggplot(data.frame(x=0,y=0))+
    geom_abline(slope = 0,intercept = 0,linetype = 2,colour = I("black"))+
    geom_step(data = fwd2,aes(x=coord,y=counts),colour = "red",direction = "vh")+
    geom_step(data = bwd2,aes(x=coord,y=counts),colour = "blue",direction = "vh")+
    theme(legend.position = "none")+
    scale_y_continuous( limits = 1.2*M * c(-1,1))+ylab(ifelse(depth > 1,"Normalized counts","Counts"))+
    scale_x_continuous(limits = c(start , end))+
    ggtitle(paste(chr,":",prettyNum(start,big.mark = ","),"-",prettyNum(end,big.mark = ",")))

  return(p)
}

stats <- stats[seqnames != "chrM"]

## library(GGally)
## stats[,log_depth := log10(depth)]
## stats[,log_width := log10(width)]


library(gridExtra)
library(RColorBrewer)
library(hexbin)

enrich <- function(data,main)
{
 
  rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
  r <- rf(16)
  
  p1 <- ggplot(data,aes(npos,depth))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+theme(legend.position = "top")+
    ylim(0,3e3)+xlim(0,500)+xlab("number of unique 5' positions mapped per region")+
    ylab("nr of reads mapped per region")

  p2 <- ggplot(data,aes(dw_ratio,pbc))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+theme(legend.position = "top")+
    xlim(0,6)+xlab("average read coverage")+ylab("unique read coverage rate")+ylim(0,1)
  
  p3 <- ggplot(data,aes(npos,width))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+theme(legend.position = "top")+
    xlim(0,500)+ylim(0,1250)+xlab("number of unique 5' positions mapped per region")+
    ylab("width of region")

    

  p4 <- ggplot(data,aes(depth,width))+stat_binhex(bins = 70)+
    scale_fill_gradientn(colours = r,trans = 'log10')+theme(legend.position = "top")+
    xlim(0,3e3)+ylim(0,1250)+ylab("width of region")+
    xlab("nr of reads mapped per region")

  grid.arrange(p1,p2,p3,p4,nrow = 2, ncol = 2,top = main)
  
}

pdf(file = file.path(figs_dir,"FoxA1-rep3_enrichment_fsr_filter.pdf"),width = 12, height = 13)
enrich(stats,"all")
enrich(stats[between(prob,.25,.75)],"fsr between .25 and .75")
enrich(stats[!between(prob,.25,.75)],"fst <  .25 or > .75")
enrich(stats[between(prob,.45,.55)],"fst between .45 and .55")
dev.off()


pdf(file = file.path(figs_dir,"FoxA1-rep3_enrichment_npos_filter.pdf"),width = 12, height = 13)
enrich(stats,"all")
enrich(stats[npos <= 5],"number of unique position <= 5")
enrich(stats[npos > 1],"number of unique position > 1")
enrich(stats[npos > 5],"number of unique position > 5")
enrich(stats[npos >10],"number of unique position > 10")
enrich(stats[npos > 25],"number of unique position > 25")
enrich(stats[npos > 50],"number of unique position > 50")
enrich(stats[npos > 100],"number of unique position > 100")
enrich(stats[npos > 200],"number of unique position > 200")
dev.off()

pdf(file = file.path(figs_dir,"FoxA1-rep3_enrichment_width_filter.pdf"),width =12,height = 13)
enrich(stats,"all")
enrich(stats[width < 100],"width <100")
enrich(stats[width > 50],"width > 50")
enrich(stats[width > 75],"width > 75")
enrich(stats[width > 100],"width > 100")
enrich(stats[width > 150],"width > 150")
enrich(stats[width > 200],"width > 200")
enrich(stats[width > 300],"width > 300")
enrich(stats[width > 500],"width > 500")
enrich(stats[width > 750],"width > 750")
dev.off()

pdf(file = file.path(figs_dir,"H3k27ac_avg_read_cover_greater_3_sample500.pdf"))
keys <- stats[dw_ratio > 3][,sample(match,500)]
lapply(keys,get_profile,stats,reads_table,fl = 1,depth = depth)
dev.off()


pdf(file = file.path(figs_dir,"H3K27ac_width_greater_500.pdf"))
keys <- stats[width > 500][,(match)]
lapply(keys,get_profile,stats,reads_table,fl = 1,depth = depth)
dev.off()


pdf(file = file.path(figs_dir,"H3K27ac_width_between_300_500_sample500.pdf"))
keys <- stats[between(width,300,500)][,sample(match,500)]
lapply(keys,get_profile,stats,reads_table,fl = 1,depth = depth)
dev.off()

pdf(file = file.path(figs_dir,"H3K27ac_width_between_100_300_sample1000.pdf"))
keys <- stats[between(width,100,300)][,sample(match,1000)]
lapply(keys,get_profile,stats,reads_table,fl = 1,depth = depth)
dev.off()

pdf(file = file.path(figs_dir,"H3K27ac_width_between_100_300_sample1000_2.pdf"))
keys <- stats[between(width,100,300)][,sample(match,1000)]
lapply(keys,get_profile,stats,reads_table,fl = 1,depth = depth)
dev.off()


## pdf(file = file.path(figs_dir,"CTCF_width_greater_300_sample500_fl50.pdf"))
## keys <- stats[width > 500][,sample(match,500)]
## lapply(keys,get_profile,stats,reads_table,fl = 50,depth = depth)
## dev.off()

pdf(file = file.path(figs_dir,"H3K27ac_width_less_100_rate_greater_50_sample500.pdf"))
keys <- stats[width < 100 & pbc > .5][,sample(match,500)]
lapply(keys,get_profile,stats,reads_table,fl = 1,depth = depth)
dev.off()


pdf(file = file.path(figs_dir,"H3K27ac_width_less_100_rate_between_25_50_sample1000_1.pdf"))
keys <- stats[width < 100 & between(pbc,.25,.5)][,sample(match,1000)]
lapply(keys,get_profile,stats,reads_table,fl = 1,depth = depth)
dev.off()

pdf(file = file.path(figs_dir,"CTCF_width_less_100_rate_less_25_dw_ratio_greater_1.pdf"))
keys <- stats[width < 100 & pbc < .25 & dw_ratio > 1][,(match)]
lapply(keys,get_profile,stats,reads_table,fl = 1,depth = depth)
dev.off()

