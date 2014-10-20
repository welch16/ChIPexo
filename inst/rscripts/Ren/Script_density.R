
rm(list = ls())
source("R/density_functions.R")
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)


datadir = "data"
load(file = file.path(datadir,"Ren_reads.RData"))


binsizes = c(300,500,750)
bins = lapply(binsizes,create.bins,genomeLength=seqlengths(reads1[[1]]))
names(bins) = as.character(binsizes)

reads = c(reads1,reads2,reads3)

counts_strand <- function(bins,reads)
{
  counts_F = countOverlaps(bins,subset(reads,subset = strand(reads) == "+"))
  counts_R = countOverlaps(bins,subset(reads,subset = strand(reads) == "-"))
  return(DataFrame(counts_F,counts_R))
}

strand_counts = mclapply(reads,function(x,bins)lapply(bins,counts_strand,x),bins,mc.cores=4)

fwd_strand_ratio <- function(f,r)return((f+1)/(f+r+2))
fwd_ratio = mclapply(unlist(strand_counts),function(x)fwd_strand_ratio(x$counts_F,x$counts_R),mc.cores=8)

ren_density = mclapply(fwd_ratio,density,mc.cores =8)


separate_name <- function(nm)
{
  split = strsplit(nm,".sort.bam.",fixed=TRUE)[[1]]
  edsn = split[1]
  bin = split[2]
  return(c(edsn,bin))
}

assign_name <- function(dens_df,name)
{
  name =  separate_name(name)
  dens_df$name = name[1]
  dens_df$bin = name[2]
  return(dens_df)
}


ren_density = mclapply(ren_density,function(x)DataFrame(x$x,x$y),mc.cores=8)

ren_density = mcmapply(assign_name,ren_density,names(ren_density),mc.cores = 8)

save(list = "ren_density",file = file.path(datadir,"RenData_density.RData"))

ren_df = do.call(rbind,mclapply(ren_density,as.data.frame,mc.cores=8))

plot_density <- function(df)
{
  p1 = ggplot(df,aes(x.x,x.y,colour = bin))+geom_line()+geom_vline(xintercept = .5,linetype="dashed")+
    facet_wrap(~name,nrow =2,scales = "free")+scale_color_brewer(palette = "Dark2")+
      theme(legend.position ="bottom")+scale_x_continuous(limits = c(0,1))+xlab("")+ylab("")
  return(p1)
}


figsdir = "inst/figs/densities"
pdf(file = file.path(figsdir,"Ren_Densities.pdf"),width = 8,height = 6)
p1 = plot_density(ren_df)
print(p1)
dev.off()
