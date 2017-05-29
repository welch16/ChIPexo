rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)
library(ggplot2)
library(RColorBrewer)

data_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files <- list.files(data_dir,pattern = "sort.bam",full.names = TRUE,include.dirs = TRUE)
files <- files[grep("bai",files,invert = TRUE)]

reads <- mclapply(files,readGAlignments,param = NULL,mc.cores = 3)
names(reads) <- gsub(".sort.bam","",basename(files))

reads <- mclapply(reads,as,"GRanges",mc.cores = 3 )

exo_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo"

fimo_files <- list.files(exo_dir,pattern = "txt",recursive = TRUE,
                         include.dirs = TRUE,full.names = TRUE)

fimo <- mclapply(fimo_files,read.table,mc.cores = 12)
fimo <- lapply(fimo,data.table)
names(fimo) <- basename(dirname(fimo_files))

##     both only_reg only_peak  rep
## 1: 12726      991       502 rep1
## 2:  3271      332       442 rep2
## 3:  4074       32      3049 rep3

# V1 motif id
# V2 sequence id
# V3 start
# V4 end
# V5 strand
# V6 score
# V7 p.value
# V8 q.value
# V9 sequence

fimo <- lapply(fimo,function(x){
  setnames(x,names(x),c("motifID","sequenceID",
                        "motifStart","motifEnd","strand",
                        "score","pval","qval","sequence"))
  return(x)})

fimo <- fimo[grep("FOXA1",names(fimo))]

chr_from_fimo <- function(x)sapply(strsplit(as.character(x),":"),
                function(y)y[1])

start_from_fimo <- function(x){
  out <- sapply(strsplit(as.character(x),":"),function(y)y[2])
  out <- sapply(strsplit(out,"-",fixed = TRUE),function(z)z[1])
  return(as.numeric(out))
}

end_from_fimo <- function(x){
  out <- sapply(strsplit(as.character(x),":"),function(y)y[2])
  out <- sapply(strsplit(out,"-",fixed = TRUE),function(z)z[2])
  return(as.numeric(out))
}

fimo <- lapply(fimo,function(x){
  x[,seqnames := chr_from_fimo(sequenceID)]
  x[,start := start_from_fimo(sequenceID)]
  x[,end := end_from_fimo(sequenceID)]
  return(x)})

fimo2 <- lapply(fimo,function(x)x[qval <= .05])

window_length <- 20
sm <- 1

## want to see if there is a spatial pre-disposition..
## we are going to divide the motifs by the strand and build a + / - 50 bp
## profile for each strand

fwd_cover <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "+")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 3)

bwd_cover <- mclapply(reads,function(x,smooth = 1){
  ff <- subset(x,as.character(strand(x)) == "-")
  ff <- resize(ff,smooth)
  return(coverage(ff))},smooth = sm,mc.cores = 3)

fwd_regions <- lapply(fimo2,function(x,window_length){
  y <- copy(x)[strand == "+"]
  out <- y[,GRanges(seqnames = seqnames,
                    ranges = IRanges(start = start + motifStart,
                      width = 1))]
  out <- resize(out,width =2*window_length + 1,fix = "center")
  return(out)
},window_length)

bwd_regions <- lapply(fimo2,function(x,window_length){
  y <- copy(x)[strand == "-"]
  out <- y[,GRanges(seqnames = seqnames,
                    ranges = IRanges(start = start + motifEnd,
                      width = 1))]
  out <- resize(out,width =2*window_length,fix = "center")
  return(out)
},window_length)


fwd_all <- mapply(function(cover,region,wl){
  mat <- cover[region]
  mat <- lapply(mat,as.vector)
  nms <- paste0(as.character(seqnames(region)),":",
                start(region),"-",end(region))
  DT <- mcmapply(function(x,nm,wl){
    data.table(coord = -wl : wl, counts = x , name = nm)
  },mat,nms,MoreArgs = list(wl),SIMPLIFY = FALSE,mc.cores = 10)
  return(do.call(rbind,DT))
},fwd_cover,fwd_regions,MoreArgs = list(window_length),SIMPLIFY = FALSE)

bwd_all <- mapply(function(cover,region,wl){
  mat <- cover[region]
  mat <- lapply(mat,as.vector)
  nms <- paste0(as.character(seqnames(region)),":",
                start(region),"-",end(region))
  DT <- mcmapply(function(x,nm,wl){
    data.table(coord = -wl : wl, counts = x , name = nm)
  },mat,nms,MoreArgs = list(wl),SIMPLIFY = FALSE,mc.cores = 10)
  return(do.call(rbind,DT))
},bwd_cover,bwd_regions,MoreArgs = list(window_length),SIMPLIFY = FALSE)


fwd_DT <- lapply(fwd_all,function(x)x[,mean(counts),by = coord])
bwd_DT <- lapply(bwd_all,function(x)x[,mean(counts),by = coord])

profile <- function(fwd,bwd,depth,repl){
  fwd <- fwd[,strand := "+"]
  bwd <- bwd[,strand := "-"]
  dt <- rbind(fwd,bwd)
##  dt[,V1 := 1e9 * V1 / depth ]
  dt[,rep :=repl]
  return(dt)
}

DT <- do.call(rbind,mapply(profile,fwd_DT,
                           bwd_DT,sapply(reads,length),
                           c("Rep-3","Rep-1","Rep-2"),SIMPLIFY = FALSE))
r <- brewer.pal(name = "Set1",n = 3)[1:2]
figs_dir <- "figs/for_paper"

pdf(file = file.path(figs_dir,"FoxA1_profiles_around_motif.pdf"),
    width= 8,height =5)
ggplot(DT,aes(coord,V1,colour = strand,linetype = rep))+
  geom_line()+
    scale_color_manual(values = rev(r))+theme_bw()+
    theme(legend.position = "top")+
  xlab("Position around motif start")+ylab("Average counts")
dev.off()

DT2 = DT
DT2[, class := paste(rep,strand,sep = "_")]
r2 <- brewer.pal(name = "Set1",n = 6)


# red #4d0000 #ff0000 #ffb3b3
# blue #00004d  #0000ff  #b3b3f
r3 = c("#00004d","#4d0000","#0000ff","#ff0000","#b3b3f","#ffb3b3")
r4 = c("darkblue","firebrick3","blue","red","lightblue","lightpink")

# dark rep1
# med rep2 
# light rep3

pdf(file = file.path(figs_dir,"FoxA1_profiles_around_motif_2.pdf"),
    width= 8,height =5)
ggplot(DT2,aes(coord,V1,colour = class))+
  geom_line(size = 1.1)+
    scale_color_manual(values = r4)+theme_bw()+
    theme(legend.position = "none")+
  xlab("Position around motif start")+ylab("Average counts")
dev.off()



## brief notes about plot

## 1 - the signal is roughly the same for rep1 and rep3
## 2 - it is lower for rep2
## 3 - these are high quality regions which binds to the motifs. there is a slight imbalance in all three replicates


## sequence plots

sequences <- lapply(fimo,function(x)x[,.(sequenceID,pval,qval,sequence)])

seqDT <- mapply(function(sequ,nm){
  motifLength <- unique(sequ[,nchar(as.character(sequence))])
  chars <- lapply(1:motifLength,function(m)sequ[,substr(sequence,m,m)])
  dt_list <- mapply(function(m,chars,sequ){
    out <- copy(sequ)
    out[,position := m]
    out[,fill := chars]
    out[,name := nm]
  },1:motifLength,chars,MoreArgs = list(sequ),SIMPLIFY = FALSE)
  return(do.call(rbind,dt_list))
},sequences,c("Rep-3","Rep-1","Rep-2"),SIMPLIFY = FALSE)

r <- brewer.pal(name = "Set1",n = 7)



pdf(file = file.path(figs_dir,"FoxA1_matched_motif_sequence.pdf"),
    width = 3,height = 6 )
mapply(function(x,y){
ggplot(x[qval < .05],aes(position,sequenceID,fill = fill))+
  geom_tile()+
  scale_fill_manual(values = r[c(3,2,6,1)],name = "")+
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(breaks = 1:11)+facet_grid( . ~ name )
},seqDT,c("Rep-3","Rep-1","Rep-2"),SIMPLIFY = FALSE)
dev.off()


outdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"

## check that the "high quality" regions overlap with peaks
peakfiles <- list.files(outdir,
   pattern = "peaks",recursive = TRUE,full.names= TRUE)
peaks <- lapply(peakfiles,
               read.table)
names(peaks) <- basename(peakfiles)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x)x[,1:3,with = FALSE])
peaks <- lapply(peaks,function(x){
  setnames(x,names(x),c("seqnames","start","end"))
  return(ChIPUtils::dt2gr(x))})

regions <- lapply(fimo,
        function(x)ChIPUtils::dt2gr(x[,.(seqnames,start,end)]))


venn <- mapply(function(region,peak){
  ov <- findOverlaps(region,peak)
  inter_peak <- length(unique(subjectHits(ov)))
  inter_reg <- length(unique(queryHits(ov)))
  out <- data.table(both = inter_reg,
                    only_reg = length(region) - inter_reg,
                    only_peak = length(peak) - inter_peak)
  return(out)                   
},regions,peaks,SIMPLIFY = FALSE)
venn <- do.call(rbind,venn)

venn[, rep := c("rep3","rep1","rep2")]

venn <- venn[order(rep)]

## only the regions were fimo found the motif
##    both only_reg only_peak  rep
## 1: 8618      529      6229 rep1
## 2: 2261      223      1915 rep2
## 3: 2742       16      4809 rep3

library(Vennerable)

pdf(file = file.path(figs_dir,"FoxA1_VennDiagram_QC_region_w_motif_and_peaks.pdf"))
venn_obj <- lapply(1:nrow(venn),function(i){
  venn[i][,Venn(SetNames = c("QC","Peaks"),
             Weight = c("11" = both,
               "10" = only_reg,"01" = only_peak))]})
lapply(venn_obj,plot)
dev.off()

### venn diagrams with all stats

load_all("~/Desktop/Docs/Code/ChIPexoQual")

data_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files <- list.files(data_dir,pattern = "sort.bam",full.names = TRUE,include.dirs = TRUE)
files <- files[grep("bai",files,invert = TRUE)]

exo <- lapply(files,create_exo_experiment,parallel = TRUE)
stats <- lapply(exo,summary_stats)

exo_reads <- lapply(exo,reads)
readlength <- mclapply(exo_reads,function(x)table(width(x)),
                mc.cores =3)

rl <- 36 ## most repeated read length

### param to get a high quality set of regions

## 1 - only regions with reads in both strands
stats <- lapply(stats,function(x)x[f > 0 & r > 0])

## 2 - wider regions
stats <- lapply(stats,function(x)x[between(width,3 * rl,1e3)])
                                   
## 3 - deeper regions
minNpos <- 15
minDepth <- 100

stats <- lapply(stats,function(x)x[npos > minNpos])
stats <- lapply(stats,function(x)x[depth > minDepth])

## remove chrM regions
stats <- lapply(stats,function(x)x[seqnames != "chrM"])

stats_gr <- lapply(stats,function(x)ChIPUtils::dt2gr(x[,2:4,with = FALSE]))
fimo_gr <- lapply(fimo,
  function(x)ChIPUtils::dt2gr(x[,.(seqnames,start,end)]))                 
                               

venn2 <- mapply(function(stat,fimo){
  ov <- findOverlaps(stat,fimo)
  inter_fimo <- length(unique(subjectHits(ov)))
  inter_stat <- length(unique(queryHits(ov)))
  out <- data.table(both = inter_stat,
                    only_stat = length(stat) - inter_stat,
                    only_fimo = length(fimo) - inter_fimo)
  return(out)                   
},stats_gr,fimo_gr,SIMPLIFY = FALSE)
venn2 <- do.call(rbind,venn2)

venn2[, rep := c("Rep-3","Rep-1","Rep-2")]
venn2 <- venn2[order(rep)]

venn2[,prop := only_stat / both ] 

##    both only_stat only_fimo   rep      prop
## 1: 7014      6614         0 Rep-1 0.9429712
## 2: 1855      1722         0 Rep-2 0.9283019
## 3: 2187      1919         0 Rep-3 0.8774577


## clearly rep1 is the best
## in its regions, more motifs had been detected

## it is questionable which is better between rep2 and rep3
## for one part the proportion of regions with motif is higher for rep2, 
## this may be adjudicated to the fact that the sample is more enriched.

## on the other hand, rep3 shows more regions with motif, even though the 
## proportion is not higher.

top_K <- function(var,K){
  return(sort(var,decreasing = TRUE)[1:min(K,length(var))])
}

all_fimo[,top_K(-log10(pval),100),by = rep]


all_fimo <- do.call(rbind,mapply(function(x,y)x[,rep := y],
  fimo,c("Rep-3","Rep-1","Rep-2"),SIMPLIFY = FALSE))                    

pdf(file = file.path(figs_dir,"FoxA1_fimo_ECDF_score.pdf"))
ggplot(all_fimo,aes(score,colour = rep))+stat_ecdf(geom = "line")+
  theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = "Replicate")+
  ylab("Empirical CDF")
ggplot(all_fimo,aes(-log10(pval),colour = rep))+stat_ecdf(geom = "line")+
  theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = "Replicate")+
  ylab("Empirical CDF")
K <- 1000
ggplot(all_fimo[,top_K(score,K),by = rep],aes(V1,colour = rep))+
  stat_ecdf(geom = "line")+xlab("score")+
  theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = "Replicate")+
  ylab("Empirical CDF")+ggtitle(paste("Top",K,"scores"))
K <- 2000
ggplot(all_fimo[,top_K(score,K),by = rep],aes(V1,colour = rep))+
  stat_ecdf(geom = "line")+xlab("score")+
  theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = "Replicate")+
  ylab("Empirical CDF")+ggtitle(paste("Top",K,"scores"))
K <- 2400
ggplot(all_fimo[,top_K(score,K),by = rep],aes(V1,colour = rep))+
  stat_ecdf(geom = "line")+xlab("score")+
  theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",name = "Replicate")+
  ylab("Empirical CDF")+ggtitle(paste("Top",K,"scores"))
dev.off()

pdf(file = file.path(figs_dir,"FoxA1_barplots.pdf"))
ggplot(venn2,aes(rep,both,fill = rep))+geom_bar(stat = "identity")+
  theme_bw()+theme(legend.position = "none")+xlab("Replicate")+
  ylab("Number of candidate sites")+
  scale_fill_brewer(palette = "Pastel1")
ggplot(venn2,aes(rep,only_stat,fill = rep))+geom_bar(stat = "identity")+
  theme_bw()+theme(legend.position = "none")+xlab("Replicate")+
  ylab("Number of sites with motif")+
  scale_fill_brewer(palette = "Pastel1")
ggplot(venn2,aes(rep,prop,fill = rep))+geom_bar(stat = "identity")+
  theme_bw()+theme(legend.position = "none")+xlab("Replicate")+
  ylab("Proportion of sites with motif")+
  scale_fill_brewer(palette = "Pastel1")+coord_cartesian(ylim = c(.75,1))
dev.off()


fimo2 <- fimo
names(fimo2) <- c("Rep-3","Rep-1","Rep-2")


topK <- c(50,250,100,500,1000,2000)

filter_topK <- function(fimo,K){
  fimo[order(-score)][1:K]
}

fimo2 <- mapply(function(x,y)x[,repl := y],fimo2,
                names(fimo2),SIMPLIFY = FALSE)

dt_list <- lapply(topK,
  function(k,fimo2){
    out <- lapply(fimo2,filter_topK,k)
    out <- do.call(rbind,out)
    out[,K := k]
    return(out)},fimo2)
dt <- do.call(rbind,dt_list)

pdf(file = "figs/FoxA1_fimo_score.pdf",width =12 ,height = 5)
ggplot(dt,aes(repl,score,fill = repl))+
  geom_boxplot(outlier.size = NA)+
  facet_grid( . ~ K)+theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 75,hjust = 1))+
  scale_fill_brewer(palette = "Pastel1")+
  xlab("Replicate")+ylab("FIMO score")
dev.off()
