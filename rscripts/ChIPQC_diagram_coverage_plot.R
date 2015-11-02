
rm(list = ls())

library(GenomicAlignments)
library(data.table)
library(ggplot2)
library(ChIPUtils)
library(scales)

file <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo/edsn931_042814_qc.sorted.bam"
param <- NULL

reads <- readGAlignments(file , param = param)

reads <- as(reads,"GRanges")

frags <- list()
frags[["fwd"]] <- reads[strand(reads) == "+"]
frags[["bwd"]] <- reads[strand(reads) == "-"]
frags[["both"]] <- reads

frag_len <- 150
frags <- lapply(frags,resize,width = frag_len)
frags <- lapply(frags,trim)


covers <- lapply(frags,coverage)

chr <- "U00096"
start <-  5000
end <- 15000

reg <- GRanges(seqnames = chr,ranges = IRanges(
  start = start, end= end),strand = "*")

covers <- lapply(covers,
  function(x)x[reg][[chr]])                

vecs <- lapply(covers,as.vector)

coord <- start:end

dt <- mapply(function(x,y)data.table( tagCounts = x, what = y),
  vecs,names(vecs),SIMPLIFY = FALSE)             

dt <- lapply(dt,function(x)x[,coords := coord])

dt <- do.call(rbind,dt)


pdf(file = "figs/for_paper/coverage_diagram.pdf",width = 12 ,height = 5)
ggplot(dt,aes(coords,tagCounts,colour = what))+geom_line()+
  scale_colour_manual(values = c("black","blue","red"))+xlab("Genomic coordinates")+
  theme_bw()+theme(legend.position = "none",axis.text = element_text(size = 0),
                   axis.ticks = element_blank())+
  geom_abline(slope = 0,intercept = 5,colour = "cornsilk4",linetype =2)+ylim(0,500)+
  ylab("Number of reads by position")
dev.off()
