

rm(list = ls())
library(viridis)
library(data.table)
library(GenomicAlignments)
library(viridis)
library(devtools)
library(parallel)
load_all("~/Desktop/Docs/Code/ChIPexoQual")

indir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files <- list.files(indir)
files <- files[grep("sort",files)]
files <- files[grep("bai",files,invert= TRUE)]
mc <- detectCores()

expts <- lapply(file.path(indir,files),create_exo_experiment,
                calc_summary = TRUE,height = 1,parallel = TRUE,mc = mc)
names(expts) <- plyr::mapvalues(gsub(".sort.bam","",files),from = c("ERR336935","ERR336942","ERR336956"),
                                to = c("Rep-3","Rep-1","Rep-2"))

stats <- lapply(expts,summary_stats)

peakdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/peaks"
peaks <- list.files(peakdir)
peaks <- lapply(file.path(peakdir,peaks),read.table)
peaks <- lapply(peaks,data.table)
names(peaks) <- names(expts)

peak_ranges <- lapply(peaks,function(x)
                      x[V8 > 200,GRanges(seqnames = V1,
                            ranges = IRanges(start = V2,end = V3))])
##stats <- lapply(stats,function(x)x[f > 0 & r > 0])
stat_ranges <- lapply(stats,function(x)ChIPUtils::dt2gr(x[,2:4,with = F]))

overlaps <- mapply(findOverlaps,peak_ranges,stat_ranges,SIMPLIFY = FALSE)

stat_peaks <- mapply(function(x,y)x[subjectHits(y)],stats,overlaps,SIMPLIFY = FALSE)
stat_peaks <- mapply(function(x,y)x[,repl := y],stat_peaks,names(stat_peaks),SIMPLIFY = FALSE)
stat_peaks <- do.call(rbind,stat_peaks)

library(scales)
library(grid)
library(gridExtra)
        

figs_dir <- "figs/for_paper"
pdf(file = file.path(figs_dir,"Carroll_mouse_FoxA1_ForwardStrandRatio_by_rep.pdf"))
ggplot(stat_peaks[f > 0 & r > 0],aes(repl,fsr,fill = repl))+geom_boxplot(outlier.shape = NA)+
  scale_fill_brewer(palette = "Pastel1")+theme_bw()+theme(legend.position = "none")+
  xlab("Replicate")+ylab("Forward Strand Ratio (FSR)")
dev.off()


### peak vs no peak

imbalance_idx <- function(p)-log(4*p*(1-p))
                         
overlap_idx <- mapply(function(stat,peak)ifelse(countOverlaps(stat,peak) > 0, "yes","no"),stat_ranges,
  peak_ranges,SIMPLIFY = FALSE)
stats <- mapply(function(x,y)x[,overlap := y],stats,overlap_idx,SIMPLIFY = FALSE)


sample_fsr <- function(stat,N = 2*nrow(stat[overlap == "yes" & f > 0 & r > 0])){
  out <- rbind(stat[overlap == "yes" & f > 0 & r > 0][sample(floor(N / 2))],
               stat[overlap == "no"][sample(floor(N / 2))])
  out
}


stats <- mapply(function(x,y)x[,repl := y],stats,names(stats),SIMPLIFY = FALSE)
stats <- do.call(rbind,stats)
setkey(stats,repl)

minD <- 30
mB <- .1
MB <- .9

## pal <- brewer_pal(3,"Pastel1")


## ggplot(stats[between(fsr,mB,MB) &  f > 0 & r > 0 & depth > minD],
##              aes(overlap , imbalance_idx(fsr),fill = repl))+
##   geom_boxplot(outlier.shape = NA)+theme_bw()+
##   theme(legend.position = "top")+
##   scale_fill_brewer(palette = "Pastel1")+
##   facet_grid( ~ repl,scales = "free_y",space = "free_y")+xlab("peak overlap")+ylab("Imbalance")+
##   coord_cartesian(ylim = c(0,.5))
## dev.off()


minD <- seq(10,200,by = 10)

alltests <- mclapply(minD,function(d,stats){
  tests <- list()
  tests[["Rep-1"]] <- stats["Rep-1"][between(fsr,mB,MB) & f > 0 & r > 0 & depth > d, wilcox.test(-log(imbalance_idx(fsr)) ~ overlap)]
  tests[["Rep-2"]] <- stats["Rep-2"][between(fsr,mB,MB) & f > 0 & r > 0 & depth > d, wilcox.test(-log(imbalance_idx(fsr)) ~ overlap)]
  tests[["Rep-3"]] <- stats["Rep-3"][between(fsr,mB,MB) & f > 0 & r > 0 & depth > d, wilcox.test(-log(imbalance_idx(fsr)) ~ overlap)]
  return(tests)},stats,mc.cores = 24)

pvals <- lapply(alltests,function(x)sapply(x,function(y)y$p.value))
pvals <- data.table(do.call(rbind,pvals))

pvals[, depth := minD]

library(reshape2)
pvals <- melt(pvals,id.vars = "depth",variable.name = "Replicate",value.name = "p.value")

save(pvals,file = "wilcoxon_pvals.RData")

pdf(file = "figs/for_paper/Carroll_FSR_depth_VS_pvalWilcoxon.pdf")
ggplot(pvals,aes(depth,-log10(p.value),colour = Replicate))+
  geom_point(size = 5)+
  geom_line(size = .5,linetype = 2)+
  xlim(50,150)+ylim(0,75)+
  theme_bw()+
  theme(legend.position = "top",
    plot.title = element_text(hjust = 0,size = 24),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22))+
  scale_color_brewer(palette = "Set1",name = "")+ggtitle("C")+
  xlab("Min. number of reads")
dev.off()
