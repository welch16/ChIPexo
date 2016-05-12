
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(gridExtra)

dr <- "/p/keles/ChIPexo/volume6/K12/downstream"
read_dr <- "/p/keles/ChIPexo/volume7/Landick/K12"
fdr <- "FDR5"
frag_len <- 50


files <- list.files(dr,recursive = TRUE,include.dirs = FALSE,
                    full.names = TRUE)
files <- files[grep(fdr,files)]

exo_files <- files[!grepl("SET",files) & !grepl("PET",files)]
exo_files <- exo_files[!grepl("motif",exo_files)]
exo_files <- exo_files[!grepl("931",exo_files) & !grepl("933",exo_files)]

set_files <- files[grepl("SET",files)]

edsn <- data.table(exo = seq(1311,1320,by = 3),
                   set = seq(1396,1402,by = 2))

exo_peaks <- lapply(exo_files[grep("peak",exo_files)],fread)
exo_sites <- lapply(exo_files[grep("sites",exo_files)],fread)

set_peaks <- lapply(set_files[grep("peak",set_files)],fread)
set_sites <- lapply(set_files[grep("sites",set_files)],fread)

dt2gr <- function(dt){
  dt <- copy(dt[,1:3,with = FALSE])
  setnames(dt,names(dt),c("seqn","st","en"))
  dt[,GRanges(seqnames = as.character(seqn),
              ranges = IRanges(start = st,
                end = en))]}
  

peaks_with_1site <- function(peak,site){
  ov <- findOverlaps(dt2gr(peak),dt2gr(site))
  peak[unique(queryHits(ov))]

}

exo_peaks_1site <- mapply(peaks_with_1site,
                          exo_peaks,exo_sites,SIMPLIFY = FALSE)

set_peaks_1site <- mapply(peaks_with_1site,
                          set_peaks,set_sites,SIMPLIFY = FALSE)

read_files <- list.files(read_dr,full.names = TRUE,recursive = TRUE,
                         pattern= "sort.bam")
read_files <- read_files[grep("bai",read_files,invert = TRUE)]
read_files <- read_files[grep("seed",read_files,invert = TRUE)]
read_files <- read_files[grep("rif",read_files)]

read_files <- read_files[grep("PET",read_files,invert = TRUE)]
read_files <- read_files[grep("1369",read_files,invert = TRUE)]

reads <- mclapply(read_files,readGAlignments,param = NULL,mc.cores = 8)
reads <- mclapply(reads,as,"GRanges",mc.cores = 8)
names(reads) <- gsub(".sort.bam","",basename(read_files))

count_strand_reads <- function(read,peak,frag_len = 1)
{
  if(frag_len > 1){
    read <- resize(read,frag_len)
  }
  f <- countOverlaps(dt2gr(peak),subset(read,strand(read) == "+"))
  r <- countOverlaps(dt2gr(peak),subset(read,strand(read) == "-"))  
  return(data.table(f ,r , fsr = f / (f + r)))
}

exo_imbalance <- mcmapply(count_strand_reads,reads[1:4],exo_peaks_1site,
                        MoreArgs = list(frag_len),
                        SIMPLIFY = FALSE,mc.cores =4 )

set_imbalance <- mcmapply(count_strand_reads,reads[5:8],set_peaks_1site,
                        MoreArgs = list(frag_len),
                        SIMPLIFY = FALSE,mc.cores =4 )

exo_imbalance <- mapply(cbind,exo_peaks_1site,exo_imbalance,SIMPLIFY = FALSE)
set_imbalance <- mapply(cbind,set_peaks_1site,set_imbalance,SIMPLIFY = FALSE)

lapply(exo_imbalance,function(x)x[,summary(f + r)])
lapply(set_imbalance,function(x)x[,summary(f + r)])

## > lapply(exo_imbalance,function(x)x[,summary(f + r)])
## [[1]]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     2.0   140.5   292.0   682.4   756.0 10410.0 

## [[2]]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      16     278     685    1229    1537   18330 

## [[3]]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     1.0   178.0   367.5  1151.0  1037.0 24450.0 

## [[4]]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      47    1673    3976    7336    9408   92130 

## > lapply(set_imbalance,function(x)x[,summary(f + r)])
## [[1]]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    69.0   634.5  1266.0  2792.0  2937.0 36500.0 

## [[2]]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      79    1124    2530    4281    5690   37670 

## [[3]]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##      77     534    1154    2796    2927   36440 

## [[4]]
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     117    1002    2197    3954    5182   34960 

qq <- .3
exo_imbalance <- lapply(exo_imbalance,
           function(x)x[f + r > quantile(f + r,prob = qq)])
set_imbalance <- lapply(set_imbalance,
           function(x)x[f + r > quantile(f + r,prob = qq)])

imbalance_dt <- function(exo,set){
  ov <- findOverlaps(dt2gr(exo),dt2gr(set))
  out <- rbind(exo[queryHits(ov),.(fsr)][,seq := "ChIP-exo"],
               set[subjectHits(ov),.(fsr)][,seq := "SE ChIP-Seq"])
  return(out)
}

dt <- mapply(imbalance_dt,exo_imbalance,set_imbalance,SIMPLIFY = FALSE)
dt <- mapply(function(x,y)x[,rif := y],dt,rep(c("0 min","20 min"),2),SIMPLIFY = FALSE)
dt <- mapply(function(x,y)x[,repl := y],dt,rep(c("Rep-1","Rep-2"),each = 2),SIMPLIFY = FALSE)
dt <- do.call(rbind,dt)

r <- brewer.pal(3,name = "Set1")

figs_dir <- "figs/supplement"

pdf(file = file.path(figs_dir,"Sig70_strand_imbalance_rif.pdf"),
    width = 6 ,height = 3 )
p1 <- ggplot(dt[grep("exo",seq)],
       aes(fsr,colour = seq))+stat_density(geom = "line")+
  facet_grid(rif ~ repl,scales = "free",space = "free")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 20,size = 6))+
  scale_color_manual(values = r[c(1)],name = "")+
  scale_x_continuous(breaks = seq(0,1,by = .25))+ylim(0,10)+
  ylab("Density")+
  xlab("Ratio of forward strand reads")
p2 <- ggplot(dt[grep("exo",seq,invert = TRUE)],
       aes(fsr,colour = seq))+stat_density(geom = "line")+
  facet_grid(rif ~ repl,scales = "free",space = "free")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 20,size = 6))+
  scale_color_manual(values = r[c(3)],name = "")+
  scale_x_continuous(breaks = seq(0,1,by = .2),limits = c(0,1))+
  ylab("Density")+
  xlab("Ratio of forward strand reads")
grid.arrange(p1,p2,nrow = 1)
dev.off()




