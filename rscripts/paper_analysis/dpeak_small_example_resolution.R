
rm(list = ls())

library(devtools)
library(parallel)
library(GenomicRanges)
library(GenomicAlignments)
library(BSgenome.Ecoli.NCBI.20080805)
library(data.table)
library(ChIPUtils)
library(reshape2)
library(RColorBrewer)

load_all("~/Desktop/Docs/Code/dpeak")

work_dir <- "/p/keles/ChIPexo/volume6/K12/other_methods/dpeak_test_mg1655"
peaks <- file.path(work_dir,"peak_to_check_dpeak.txt")
reads <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic/edsn931_Sig70.sort.bam"

mc <- detectCores()

system.time(
dpeak <- dpeakRead(peakfile = peaks,readfile = reads,
                  fileFormat = "bam",parallel = TRUE,nCore = mc,fragLen = 50)
            )

## note there is a bug from GRanges where no peak can start at zero
## there is an issue with dpeak:
## for one side

peak_tmp <- tempfile(pattern = "peak_",fileext = ".txt")
aux <- data.table(read.table(peaks))
aux[, V1 := "NC_010473"]
aux[, V1 := "NC_000913"]
write.table(aux,file = peak_tmp,quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)

memeArg <- paste("-dna -mod zoops -nmotifs 1 -minw 10 -maxw 20 -revcomp -maxsize 1000000000 -nostatus","-p",mc,sep = " ")
fimoArg <- "-max-stored-scores 100000000 -motif-pseudo 0.000001 --verbosity 1"

system.time(
motif <- dpeakMotif(peakfile = peak_tmp,refGenome = Ecoli,
  memeArgument = memeArg,fimoArgument = fimoArg,
  tempDir = work_dir)
)
            

system.time(
fit <- dpeakFit(dpeak,objectMotif = motif,nCore = mc,maxComp = 5)
)
            
##     user   system  elapsed 
## 3648.378   35.019  305.228 

gem_sites <- data.table(
  read.table("/p/keles/ChIPexo/volume6/K12/other_methods/dpeak_test_mg1655/gem_example/gem_example_outputs/gem_example_1_GEM_events.bed",skip = 1))                        

gemotif <- motif


motpeaks <- IRanges(start = motif@peakStart,end = motif@peakEnd)
motpeaks <- split(motpeaks,1:length(motpeaks))

gems <- lapply(motpeaks,function(x,gem_sites){
  ir <- IRanges(start = gem_sites$V2,end = gem_sites$V3)
  ir <- IRanges(start  = mid(ir), width = 1)
  return(sort(mid(subsetByOverlaps(ir , x))))
},gem_sites)
names(gems) <- NULL

gemotif@locMotif <- gems

system.time(
fit2 <- dpeakFit(dpeak,objectMotif = gemotif,nCore = mc,maxComp = 5)
)

save(fit,fit2,file = "dpeakFit_motif_gem_init_mg1655.RData")


reads <- readGAlignments(file = reads, param = NULL)
reads <- as(reads,"GRanges")
peaks <- data.table(read.table(peaks))
peaks <- peaks[,GRanges(seqnames = V1, ranges = IRanges(start = V2,end = V3))]
peaks <- split(peaks,1:length(peaks))

base_dir <- "/p/keles/ChIPexo/volume6/K12"

annot_dir <- file.path(base_dir,"annotations")
afiles <- list.files(annot_dir)
afiles <- afiles[grep("csv",afiles)]
annots <- read.csv(file.path(annot_dir,afiles),sep = "\t")
annots <- data.table(annots)
annots <- annots[grep("Sigma70",sigma.factor)]
annots <- annots[grep("TIM",description)]
annots <- annots[grep("RPP",description,invert = TRUE)]
annots <- annots[grep("AIPP",description,invert = TRUE)]
annots <- annots[grep("WHO",description,invert = TRUE)]
annots <- annots[grep("NTAS",description,invert = TRUE)]
annots <- annots[grep("ICA",description,invert = TRUE)]
annotDT <- annots
annots <- annots[,IRanges(start = coord,width = 1)]
annots <- sort(annots)
annots <- mclapply(peaks,function(x)mid(subsetByOverlaps(annots,ranges(x))),mc.cores = mc)
names(annots) <- NULL

plot_peak <- function(peak,reads,fl = 1)
{
  fwd <- reads[strand(reads) == "+"]
  bwd <- reads[strand(reads) == "-"]
  width(fwd) <- fl
  width(bwd) <- fl
  fwdc <- coverage(ranges(fwd))[ranges(peak)]
  bwdc <- coverage(ranges(bwd))[ranges(peak)]
  dt <- data.table( coord = start(peak):end(peak),
                   fwd = as.vector(fwdc),
                   bwd = - as.vector(bwdc))
  dt <- melt(dt, id.vars = "coord")
  lim <- max(abs(dt[,range(value)]))*1.2

  out <- ggplot(dt , aes(coord ,value , colour = variable ))+geom_step()+
    scale_color_brewer(palette = "Set1",name = "")+theme_bw()+
    theme(legend.position = "none")+ylim(-lim,lim)+
    xlab("Genomic coordinate")+ylab("Counts")+xlim(start(peak),end(peak))

  return(out)

}

add_annot <- function(peak_cover,annot , col )
{
  out <- peak_cover+geom_vline(xintercept = annot,colour = col,linetype = 2,size = 1)
  return(out)
}

figs_dir <- "/p/keles/ChIPexo/volume3/ChIPexo/figs/gem_dpeak_reso"

peaks_coverage <- mclapply(peaks,plot_peak,reads,mc.cores = mc)

peaks_coverage <- mapply(add_annot,peaks_coverage,annots,MoreArgs = list("black"),
                         SIMPLIFY = FALSE)

colors <- brewer.pal(8,name = "Dark2")

pdf(file = file.path(figs_dir,"profile_with_annot_mg1655.pdf"),width = 9,height = 6)
u <- lapply(peaks_coverage,print)
dev.off()

pdf(file = file.path(figs_dir,"profile_motif_init_mg1655.pdf"),width = 9 , height = 6)
u <- lapply(mapply(add_annot,peaks_coverage,motif@locMotif,MoreArgs = list(colors[1]),
                         SIMPLIFY = FALSE),print)
dev.off()

pdf(file = file.path(figs_dir,"profile_gem_init_mg1655.pdf"),width = 9,height = 6)
u <- lapply(mapply(add_annot,peaks_coverage,gemotif@locMotif,MoreArgs = list(colors[2]),
                         SIMPLIFY = FALSE),print)
dev.off()

pdf(file = file.path(figs_dir,"profile_dpeak_est_motif_mg1655.pdf"),width = 9, height = 6)
u <- lapply( mapply(add_annot,peaks_coverage,fit@optMu,MoreArgs = list(colors[3]),
                         SIMPLIFY = FALSE),print)
dev.off()

pdf(file = file.path(figs_dir,"profile_dpeak_est_gem_mg1655.pdf"),width = 9,height = 6)
u <- lapply(mapply(add_annot,peaks_coverage,fit2@optMu,MoreArgs = list(colors[4]),
                         SIMPLIFY = FALSE),print)
dev.off()


library(gridExtra)

pdf(file = file.path(figs_dir,"profile_dpeak_est_comp_mg1655.pdf"),width = 9, height = 7)
u <- mapply(grid.arrange,mapply(add_annot,peaks_coverage,fit@optMu,MoreArgs = list(colors[3]),
                         SIMPLIFY = FALSE),
            mapply(add_annot,peaks_coverage,fit2@optMu,MoreArgs = list(colors[4]),
                         SIMPLIFY = FALSE),MoreArgs = list(nrow = 2))
dev.off()


pdf(file = file.path(figs_dir,"profile_motif_init_comp_mg1655.pdf"),width = 9, height = 7)
u <- mapply(grid.arrange,mapply(add_annot,peaks_coverage,motif@locMotif,MoreArgs = list(colors[1]),
                         SIMPLIFY = FALSE),
            mapply(add_annot,peaks_coverage,fit@optMu,MoreArgs = list(colors[3]),
                         SIMPLIFY = FALSE),MoreArgs = list(nrow = 2))
dev.off()

pdf(file = file.path(figs_dir,"profile_gem_init_comp_mg1655.pdf"),width = 9, height = 7)
u <- mapply(grid.arrange,mapply(add_annot,peaks_coverage,gemotif@locMotif,MoreArgs = list(colors[2]),
                         SIMPLIFY = FALSE),
            mapply(add_annot,peaks_coverage,fit2@optMu,MoreArgs = list(colors[4]),
                         SIMPLIFY = FALSE),MoreArgs = list(nrow = 2))
dev.off()


pdf(file = file.path(figs_dir,"profile_init_strat_comp_mg1655.pdf"),width = 9, height = 7)
u <- mapply(grid.arrange,mapply(add_annot,peaks_coverage,motif@locMotif,MoreArgs = list(colors[1]),
                         SIMPLIFY = FALSE),
            mapply(add_annot,peaks_coverage,gemotif@locMotif,MoreArgs = list(colors[2]),
                         SIMPLIFY = FALSE),MoreArgs = list(nrow = 2))
dev.off()


## estimated param analysis:

delta <- rbind(data.table(method = "motif",delta = do.call(c,fit@optDelta)),
                data.table(method = "gem",delta = do.call(c,fit2@optDelta)))
sigma <- rbind(data.table(method = "motif",sigma = do.call(c,fit@optSigma)),
               data.table(method = "gem",sigma = do.call(c,fit2@optSigma)))
pi0 <- rbind(data.table(method = "motif",pi0 = do.call(c,fit@optPi0)),
             data.table(method = "gem",pi0 = do.call(c,fit2@optPi0)))

pdf(file = file.path(figs_dir,"param_plot_method_mg1655.pdf"))
ggplot(delta,aes(method , delta))+geom_boxplot()
ggplot(sigma,aes(method , sigma))+geom_boxplot()
ggplot(pi0,aes(method , pi0))+geom_boxplot()+theme_bw()+ylim(0,.3)+
  theme(axis.text = element_text(size = 12),legend.position ="none")
dev.off()


## resolution of the two methods?
resolution <- function(annot,fit)
{
  out <- sort(sapply(annot,function(x)min(abs(x - fit))))
  n <- length(annot)
  m <- length(fit)
  if(n >= m){
    out <- out[1:m]
  }
  return(out)   
}

dpeak_motif_reso <- mapply(resolution,annots,fit@optMu,SIMPLIFY = FALSE)
dpeak_gem_reso <- mapply(resolution,annots,fit2@optMu,SIMPLIFY = FALSE)
motif_reso <- mapply(resolution,annots,motif@locMotif,SIMPLIFY = FALSE)
gem_reso <- mapply(resolution,annots,gemotif@locMotif,SIMPLIFY  = FALSE)

reso <- rbind(data.table(method = "dpeak (motif)",reso = do.call(c,dpeak_motif_reso)),
              data.table(method = "dpeak (gem)",reso = do.call(c,dpeak_gem_reso)),
              data.table(method = "gem",reso = do.call(c,gem_reso)))

              data.table(method = "motif",reso = do.call(c,motif_reso))

pdf(file = file.path(figs_dir,"resolution_mg1655.pdf"))
ggplot(reso,aes(method,reso))+geom_boxplot()+
  theme_bw()+theme(legend.position = "none")+ylab("Resolution")+
  ylim(0,30)
dev.off()





cc <- coverage(reads)

zero_ranges <- function(dd,minW = 1,minGap = 10)
{
  n <- length(dd)
  id <- which(dd > 0)
  b1 <- c(1,id + 1)
  b2 <- c(id - 1,n)
  w <- b2 - b1 + 1
  rr <- IRanges(start  =b1[w > minW - 1]  , end =  b2[w > minW - 1])
  return(reduce(rr,min.gapwidth = minGap)) 
}


library(RANN)



pdf(file.path(figs_dir,"profiles_chipexo_filter_regions.pdf"),width = 12, height = 4)
plots <- mclapply(1:length(peaks),function(h){
  cc_h <- cc[peaks[[h]]][[1]]
  dt_h <- data.table(as.vector(cc_h))
  A <-  nn2(dt_h,dt_h,k = 3)
  p1 <- plot_peak(peaks[[h]],reads,58)
  z <- zero_ranges(A[[2]][,2],minW = 10,minGap = 20)
  if(length(z) >= 3){
    z <- z[width(z) > quantile(width(z),prob = .5)]
  }
  out <- add_annot(p1,annots[[h]],"black")
  out <- add_annot(out,motif@locMotif[[h]],colors[1])
  if(length(z) > 0){
    z <- shift(z,start(peaks[[h]]) - 1)
    z1 <- end(z)
    if(z1 > end(peaks[[h]])){
      z1 <- end(peaks[[h]])
    }
    suppressWarnings(
    out <- out + annotate("rect",xmin = start(z) ,
             xmax = z1 ,
             ymin = -Inf , ymax = Inf,colour = "grey",fill = "orange",alpha = .25)
                     )
  }  
  return(out)}
  ,mc.cores = mc)
u <- lapply(plots,print)
dev.off()

