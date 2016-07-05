
rm(list = ls())
library(viridis)
library(data.table)
library(GenomicAlignments)
library(viridis)
library(devtools)
library(parallel)
load_all("~/Desktop/Docs/Code/ChIPexoQual")

indir <- "data/ChIPexo_QC_runs"
files <- list.files(indir,full.names = TRUE)

files <- files[grepl("carr",files) & grepl("mous",files)]

mc <- detectCores()

load_stat <- function(x)
{
  load(x)
  ext_stats$stats
}


stats <- mclapply(files,load_stat,mc.cores = mc)
names(stats) <- paste("Rep",1:3,sep = "-")
peakdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/peaks"
peaks <- list.files(peakdir)
peaks <- lapply(file.path(peakdir,peaks),fread)
names(peaks) <- names(stats)

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
ggplot(stat_peaks[f > 0 & r > 0],aes(repl,-log10(4* fsr * (1 - fsr)),fill = repl))+geom_boxplot(outlier.shape = NA)+
  scale_fill_brewer(palette = "Pastel1")+theme_bw()+theme(legend.position = "none")+
  xlab("Replicate")+ylab("Forward Strand Ratio (FSR)")+ylim(0,.2)
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

save(pvals,file = "data/wilcoxon_pvals.RData")

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

library(scales)


wilcoxon_dt <- function(minD,stats)
{
  dt <- stats[,.(depth,fsr,imb,overlap,repl)]
  dt <- dt[depth > minD]

  dt[,dd := minD]

  return(dt)
}



wilcoxon_plots <- function(minD,stats,imbalance = FALSE)
{
  dt <- stats[,.(depth,fsr,imb,overlap,repl)]
  dt <- dt[depth > minD]

  out <- ggplot(dt,aes_string("overlap",
                              ifelse(imbalance,"imb","fsr"),
                              fill = "repl"))+
      geom_boxplot(outlier.size = NA)+theme_bw()+
      theme(legend.position = "top")+
      scale_fill_brewer(palette = "Pastel1")+
      ggtitle(minD)+facet_wrap(  ~ repl,nrow = 1)
  if(imbalance){
    out <- out + scale_y_log10()+ylab("imbalance index")
  }else{
    out <- out + ylim(0,1)+ylab("forward strand ratio")
  }
  print(out)
}




stats[,imb := imbalance_idx(fsr)]


minD <- seq(50,200,by = 10)
dt_list <- mclapply(minD,wilcoxon_dt,stats,mc.cores = 20)

dt <- do.call(rbind,dt_list)

pdf(file = "figs/strand_imbalance_FoxA1.pdf",width = 12,height = 5)
ggplot(dt,aes(repl,imb,fill = overlap))+  
      geom_boxplot(outlier.size = NA,width = .8,position = "dodge")+
  theme_bw()+
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 75,hjust  = 1))+
      scale_fill_brewer(name = "Overlap",palette = "Pastel2")+
      scale_y_log10(limits = c(1e-5,1e1))+
  ylab("Strand Imbalance")+facet_grid(. ~ dd)+
  xlab("Replicate")
ggplot(dt[dd %in% c(50,100,150,200)],aes(repl,imb,fill = overlap))+  
      geom_boxplot(outlier.size = NA,width = .8,position = "dodge")+
  theme_bw()+
      theme(legend.position = "top",
            axis.text.x = element_text(angle = 75,hjust  = 1))+
      scale_fill_brewer(name = "Overlap",palette = "Pastel2")+
      scale_y_log10(limits = c(1e-5,1e1))+
  ylab("Strand Imbalance")+facet_grid(. ~ dd)+
  xlab("Replicate")
dev.off()


pdf("strand_ratio.pdf")
plots <- lapply(minD,wilcoxon_plots,stats,imbalance = FALSE)
dev.off()
pdf("imbalance.pdf")
plots <- lapply(minD,wilcoxon_plots,stats,imbalance = TRUE)
dev.off()

library(dplyr)

clean_test <- function(minD,stats,repli,lb,ub){

  dt <- stats[repli] %>% filter( f > 0 & r > 0) %>%
    filter(depth > minD) %>% filter(between(fsr,lb,ub))
  dt <- select(dt,imb,fsr,overlap)

  tests <- list()
  tests[[1]] <- dt[,broom::tidy(wilcox.test( imb ~ overlap))]
  
  imb1 <- dt %>% filter(overlap == "yes") %>% select(imb)
  imb2 <- dt %>% filter(overlap == "no") %>% select(imb)
  
  tests[[2]] <- broom::tidy( ks.test(imb1[[1]],imb2[[1]]))

  tests <- lapply(tests,data.table)

  names(tests) <- c("wilcox","ks")

  tests <- mapply(function(x,y)x[, test := y ],tests,names(tests),
                  SIMPLIFY = FALSE)
  tests <- do.call(rbind,tests)

  tests[,repl := repli]
  tests[,depth := minD]
                  
  return(tests)
  
  
}

minD <- seq(10,500,by = 10)

rep1 <- mclapply(minD,clean_test,stats,"Rep-1",.1,.9,
                 mc.cores = mc) 
rep2 <- mclapply(minD,clean_test,stats,"Rep-2",.1,.9,
                 mc.cores = mc) 
rep3 <- mclapply(minD,clean_test,stats,"Rep-3",.1,.9,
                 mc.cores = mc)

tests <- do.call(rbind,
  list(do.call(rbind,rep1),
       do.call(rbind,rep2),
       do.call(rbind,rep3)))

save(tests,file = "data/carroll_mouse_imbalance_test.RData")





