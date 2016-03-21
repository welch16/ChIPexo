
rm(list = ls())
library(GenomicAlignments)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(scales)
library(viridis)

fimodir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/fimo"
fimofiles <- list.files(fimodir,recursive = TRUE)

fimofiles <- fimofiles[grep("txt",fimofiles)]
fimofiles <- fimofiles[grep("FOXA1",fimofiles)]

fimo <- lapply(file.path(fimodir,fimofiles),function(x){
  message(x)
  read.table(x,stringsAsFactors = FALSE)})
fimo <- lapply(fimo,function(x){
  x <- data.table(x)
  setnames(x,names(x),c("pattern","sequence","start","stop","strand","score",
                        "p.value","q.value","matched"))
  return(x)})
names(fimo) <- gsub(".meme","",dirname(fimofiles),fixed = TRUE)

fimo <- mapply(function(x,y){
  z <- strsplit(y,"_",fixed = TRUE)[[1]]  
  x[, rep := z[1] ]
  x[, TF := z[2] ]
  x[, db := do.call(paste0,as.list(z[-c(1:2)]))]
  return(x)},fimo,names(fimo),SIMPLIFY = FALSE)

fimo <- do.call(rbind,fimo)

sequences <- fimo[,(sequence)]
sequences <- strsplit(sequences,":",fixed = TRUE)

fimo[ ,sitechrID := sapply(sequences,function(x)x[1])]
sequences <- strsplit(sapply(sequences,function(x)x[2]),"-",fixed = TRUE)

fimo[,sitestart := sapply(sequences ,function(x)as.numeric(x[1]))]
fimo[,siteend := sapply(sequences ,function(x)as.numeric(x[2]))]
rm(sequences)

fimo[,center := mid(IRanges(start = start, end = stop))]

fimo[,length(pattern),by = .(rep,TF,db)]

DT <- fimo[TF == "FOXA1"]
DT[, rep := plyr::mapvalues(rep,from = c("ERR336935","ERR336942","ERR336956"),
       to = c("Rep-3","Rep-1","Rep-2"))]

## get sites

sitedir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/sites"
files <- list.files(sitedir)

sites <- lapply(file.path(sitedir,files),read.table,header = TRUE)
sites <- lapply(sites,data.table)

peakdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/peaks"
peaks <- list.files(peakdir)
peaks <- lapply(file.path(peakdir,peaks),read.table)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x){
  setnames(x,names(x),c("chrID","peakStart",
                        "peakStop","peakSize",
                        "logAveP","logMinP","aveLogP",
                        "aveChipCount","maxChipCount","map","GC"))
  x})
peaks <- lapply(peaks,function(x)x[,peakID := paste0(chrID,":",peakStart,"-",peakStop)])

lapply(peaks,function(x)x[,summary(aveChipCount)])

peaks <- lapply(peaks,function(x)x[aveChipCount > 250])

ids <- lapply(peaks,function(x)x[,(peakID)])

sites <- mapply(function(x,y)x[peakID %in% y],sites,ids,SIMPLIFY= FALSE)
names(sites) <- sapply(strsplit(files,"_",fixed = TRUE),function(x)x[1])

sites <- mapply(function(x,y)x[,rep := y],sites,names(sites),SIMPLIFY = FALSE)
sites <- do.call(rbind,sites)
sites[, rep := plyr::mapvalues(rep,from = c("ERR336935","ERR336942","ERR336956"),
       to = c("Rep-3","Rep-1","Rep-2"))]



##

DTsummary <- DT[TF == "FOXA1",.(nMotif = length(pattern)),by = rep]
DTsummary <- merge(DTsummary,sites[,.(nSites = length(strength)), by = rep],by = "rep")
DTsummary[, Motif_perc := nMotif / nSites]

### proportion plot

figs_dir <- "figs/FOXA1_mm9_fimo"
dir.create(figs_dir,showWarnings = FALSE)

pdf(file = file.path(figs_dir,"FOXA1_summaries_by_replicate.pdf"))
ggplot(DTsummary,aes(rep,nMotif,fill = rep))+
  geom_bar(stat = "identity",colour = "black")+
  theme_bw()+
  theme(legend.position = "none",plot.title = element_text(hjust = 0,size = 32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26))+
  xlab("Replicates")+ylab("Number of Motif occurrences")+
  scale_fill_brewer(palette = "Pastel1")+ggtitle("B")
ggplot(DTsummary,aes(rep,nSites,fill = rep))+
  geom_bar(stat = "identity",colour = "black")+
  theme_bw()+
  theme(legend.position = "none",plot.title = element_text(hjust = 0,size = 32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26))+
  xlab("Replicates")+ylab("Number of FoxA1 binding events")+
  scale_fill_brewer(palette = "Pastel1")
ggplot(DTsummary,aes(rep,Motif_perc,fill = rep))+
  geom_bar(stat = "identity",colour = "black")+
  scale_y_continuous(labels = percent)+
  theme_bw()+
  theme(legend.position = "none",plot.title = element_text(hjust = 0,size = 32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26))+  
  xlab("Replicates")+ylab("Percentage of sites with FoxA1 motif")+
  scale_fill_brewer(palette = "Pastel1")+ggtitle("C")
dev.off()


### cummulative plots

pdf(file = file.path(figs_dir,"FOXA1_log10pval_density.pdf"))
ggplot(DT,aes(-log10(p.value),fill = rep))+
  geom_histogram(aes(y = ..density..),colour = "black")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_fill_brewer(palette = "Pastel1")+
  facet_grid( rep ~ . ,scales = "free_y")+
  scale_x_continuous(limits = c(4,5.5))+ylab("Density")
dev.off()


## x-axis =  -log10(p.value)
## y-axis = nr (or prop) of motifs with p.value below that

pvals <- sort(DT[,(p.value)])
nmotifs <- lapply(pvals,function(p,DT,DTsummary){
  base <- data.table(rep = paste0("Rep-",1:3), cumMotif = 0)
  out <- DT[p.value <=  p][,.(cumMotif = length(pattern)),by = rep]
  if(nrow(out) == 0){
    out <- base
  }else if(nrow(out) < 3){
    setkey(base,rep)
    suppressWarnings(base[ out[,(rep)] , cumMotif := out[,(cumMotif)] ])
  }
  out[, p := p]
  out <- merge(out,DTsummary,by = "rep")
  out[, motif_perc_cumm:= cumMotif / nSites]
  out[,det_motif_perc_cumm := cumMotif / nMotif]
  out},DT,DTsummary)
nmotifs <- do.call(rbind,nmotifs)


pdf(file = file.path(figs_dir,"FOXA1_log10pval_cumulative_nmotifs.pdf"))
p <- ggplot(nmotifs,aes(-log10(p) , motif_perc_cumm,colour = rep))+
  geom_line(size = 1)+theme_bw()+
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        plot.title = element_text(hjust = 0,size = 32),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26),
        legend.text = element_text(size = 20))+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = ""))+
  scale_y_continuous(labels = percent)+
  xlab("-log10(p.value)")+
  ylab("Percentage of detected motifs")+ggtitle("D")
p
p + scale_x_continuous(limits = c(4,5))
dev.off()

pdf(file = file.path(figs_dir,"FOXA1_log10pval_cumulative_nmotifs_detect.pdf"))
p <- ggplot(nmotifs,aes(p , det_motif_perc_cumm,colour = rep))+geom_line()+theme_bw()+
  theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1",guide = guide_legend(title = "Replicate"))+
  scale_y_continuous(labels = percent)+xlab("-log10(p.value)")+
  ylab("Percentage of detected motifs")
p
p + scale_x_continuous(limits = c(4,5))
dev.off()

                  
## 1 - Rep-2 is quite bad
## 2 - Rep-1 and Rep-3 are equivalente in terms of statistical significance of motifs found
##     but this analysis is restringed to the fact that we are finding more motifs for Rep-1

## x-axis = distance in bp
## y-axis = nr (or prop) of binding events within that distance

## need to calculate distances
