
rm(list = ls())
source("R/density_functions.R")
library(ggplot2)
library(GenomicAlignments)
library(parallel)
library(knitr)

load("data/chip.exo.RData")
load("data/chip.seq.pet.RData")
load("data/sample.summary.RData")
load("data/chip.seq.set.from.pet.RData")

binsizes = c(300,500,750)
bins = lapply(binsizes,create.bins,genomeLength=seqlengths(exo[[1]]))
names(bins) = as.character(binsizes)

counts_strand <- function(bins,reads)
{
  counts_F = countOverlaps(bins,subset(reads,subset = strand(reads) == "+"))
  counts_R = countOverlaps(bins,subset(reads,subset = strand(reads) == "-"))
  return(DataFrame(counts_F,counts_R))
}

exo_strand_counts = mclapply(exo,function(x,bins)lapply(bins,counts_strand,x),bins,mc.cores=6)
pet_strand_counts = mclapply(pet,function(x,bins)lapply(bins,counts_strand,x),bins,mc.cores=6)
set_strand_counts = mclapply(set,function(x,bins)lapply(bins,counts_strand,x),bins,mc.cores=6)

fwd_strand_ratio <- function(f,r)return((f+1)/(f+r+2))
exo_fwd_ratio = mclapply(unlist(exo_strand_counts),function(x)fwd_strand_ratio(x$counts_F,x$counts_R),mc.cores=8)
pet_fwd_ratio = mclapply(unlist(pet_strand_counts),function(x)fwd_strand_ratio(x$counts_F,x$counts_R),mc.cores=8)
set_fwd_ratio = mclapply(unlist(set_strand_counts),function(x)fwd_strand_ratio(x$counts_F,x$counts_R),mc.cores=8)


exo_density = mclapply(exo_fwd_ratio,density,mc.cores =8)
pet_density = mclapply(pet_fwd_ratio,density,mc.cores =8)
set_density = mclapply(set_fwd_ratio,density,mc.cores =8)

exo_density = mclapply(exo_density,function(x)DataFrame(x$x,x$y),mc.cores=8)
pet_density = mclapply(pet_density,function(x)DataFrame(x$x,x$y),mc.cores=8)
set_density = mclapply(set_density,function(x)DataFrame(x$x,x$y),mc.cores=8)

separate_name <- function(nm)
{
  split = strsplit(nm,"_")[[1]]
  edsn = split[1]
  split = strsplit(split[2],".",fixed=TRUE)[[1]]
  bin = split[2]
  return(c(edsn,bin))
}

assign_name <- function(dens_df,name)
{ 
  split = separate_name(name)
  dens_df$edsn = split[1]
  dens_df$bin = split[2]
  return(dens_df)
}

exo_density = mcmapply(assign_name,exo_density,names(exo_density),mc.cores = 8)
pet_density = mcmapply(assign_name,pet_density,names(pet_density),mc.cores = 8)
set_density = mcmapply(assign_name,set_density,names(set_density),mc.cores = 8)


save(list = c("exo_density","pet_density","set_density"),file = "data/densities.RData")


exo_density = do.call(rbind,mclapply(exo_density,as.data.frame,mc.cores=8))
pet_density = do.call(rbind,mclapply(pet_density,as.data.frame,mc.cores=8))
set_density = do.call(rbind,mclapply(set_density,as.data.frame,mc.cores=8))

plot_density <- function(df)
{
  p1 = ggplot(df,aes(x.x,x.y,colour = bin))+geom_line()+geom_vline(xintercept = .5,linetype="dashed")+
    facet_wrap(~edsn,ncol=4,scales = "free")+scale_color_brewer(palette = "Dark2")+
      theme(legend.position ="bottom")+scale_x_continuous(limits = c(0,1))+xlab("")+ylab("")
  return(p1)
}

exo_density_1 = exo_density[grepl("13",exo_density$edsn),]
exo_density_2 = exo_density[!grepl("13",exo_density$edsn),]


p1 = plot_density(exo_density_1)
p2 = plot_density(exo_density_2)
p3 = plot_density(pet_density)
p4 = plot_density(set_density)

figsdir = "inst/figs/densities"
pdf(file = file.path(figsdir,"Densities.pdf"),width = 8,height = 12)
print(p1 + ggtitle("ChIP-exo 13-- fwd. strand ratio densities"))
print(p2 + ggtitle("ChIP-exo 9-- fwd. strand ratio densities"))
print(p3 + ggtitle("ChIP-seq PET fwd. strand ratio densities"))
print(p4 + ggtitle("ChIP-seq SET fwd. strand ratio densities"))
dev.off()

