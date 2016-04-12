#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    arc_vs_urcr_factor_test.R - Generates the arc vs urcr plots with total depth, width and npos factors
      in both statistics

  Arguments:

   -- infile

      RData file with the QC pipeline statistics

   -- outdir

      Dri where the pdf file with plots is gonna be saved

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 2)

library(parallel)
library(data.table)
library(ggplot2)
library(scales)
library(hexbin)
library(viridis)
library(gridExtra)

infile <- args[1]
outdir <- args[2]

load(infile)

stats <- ext_stats[[1]]

single_stats <- stats[f > 0 & r > 0 ]

get_factors <- function(stat){
  all_depth <- stat[,sum(depth)]
  all_width <- stat[,sum(width)]
  all_npos <- stat[,sum(npos)]

  out <- c(all_depth / all_width,
           all_depth / all_npos,
           all_npos / all_width)
  return(out)
}

factors <- get_factors(stats)
single_factors <- get_factors(single_stats)

r <- viridis(1e3, option = "D")

stats <- stats[,.(ave_reads, cover_rate)]
single_stats <- single_stats[,.(ave_reads,cover_rate)]

stats1 <- copy(stats)[,ave_reads := ave_reads / factors[1] ] 
single_stats1 <- copy(single_stats)[,ave_reads := ave_reads / single_factors[1] ]

stats2 <- copy(stats)[,ave_reads := ave_reads / factors[2] ] 
single_stats2 <- copy(single_stats)[,ave_reads := ave_reads / single_factors[2] ]

stats3 <- copy(stats)[,ave_reads := ave_reads * factors[3] ] 
single_stats3 <- copy(single_stats)[,ave_reads := ave_reads * single_factors[3] ]

outfile <- file.path(outdir,gsub(".RData",".pdf",basename(infile)))

pdf(file = outfile,width = 6,height = 6)
ggplot(stats,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")+xlim(0,12)+
  ggtitle("All, current")
ggplot(stats1,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")+xlim(0,20)+
  ggtitle("All, arc = depth / DEPTH * WIDTH / width")
ggplot(stats2,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")+xlim(0,8)+
  ggtitle("All, arc = depth / DEPTH * NPOS / width")
ggplot(stats3,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")+xlim(0,1)+
  ggtitle("All, arc = depth / WIDTH * NPOS / width")
ggplot(single_stats,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")+xlim(0,12)+
  ggtitle("Single, current")
ggplot(single_stats1,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")+xlim(0,20)+
  ggtitle("Single, arc = depth / DEPTH * WIDTH / width")
ggplot(single_stats2,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")+xlim(0,8)+
  ggtitle("Single, arc = depth / DEPTH * NPOS / width")
ggplot(single_stats3,aes(ave_reads,cover_rate))+stat_binhex(bins = 70)+
  theme_bw()+theme(legend.position = "top")+
  scale_fill_gradientn(colours = r,trans = 'log10',labels=trans_format('log10',math_format(10^.x)) )+  
  xlab("Average read coverage (ARC)")+ylab("Unique read coverage rate (URCR)")+xlim(0,1)+
  ggtitle("Single, arc = depth / WIDTH * NPOS / width")
dev.off()




