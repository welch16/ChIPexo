 rm( list=ls() )

graphics.off()
library(dpeak)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(mosaics)

# options

extension = 20

#dir.peak <- "/scratch/deconvolution/comparison_qSET/"
dir.out = "figs/for_paper/"

# read promoter regions for sigma70 from regulonDB

prom = read.table(
    "/p/keles/ENCODE2Data/volume4/KelesGroup_DongjunChung/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/database_regulonDB/PromoterSigma70Set.txt",
  skip=33, sep="\t")
loc.prom = prom[,4]
prom = data.table(prom)

# load peak region & motif
load("/p/keles/ChIPexo/volume6/K12/subsample_resolution_sensitivity/subsample_reso_peaks.RData") ## peaks

peaks = unlist(peaks)

pfiles = vapply(seq_along(peaks),function(x)tempfile(fileext = ".bed"),"")

read_peak <- function(peak,filename)
{
  export(peak,type = "bed",filename = filename)
  dt =  fread(filename)
  dt[,V2 := V2 + 1]
  dt[,V3 := V3 + 1]
  dt
}

exoidx = 1:4
petidx = 5:6
setidx = 7:8

peaks_bed = mapply(read_peak,peaks,pfiles,SIMPLIFY = FALSE)
names(peaks_bed) = names(peaks)

peakList.exo = peaks_bed[exoidx]
peakList.PET = peaks_bed[petidx]
peakList.SET = peaks_bed[setidx]

# load predictions

load("/p/keles/ChIPexo/volume6/K12/subsample_resolution_sensitivity/ChIPexo/sites/TFBS_separate.RData")
fit.exo_list = exofit

load("/p/keles/ChIPexo/volume6/K12/subsample_resolution_sensitivity/ChIPseq_PET/sites/TFBS_separate.RData")
fit.PET_list = petfit

load("/p/keles/ChIPexo/volume6/K12/subsample_resolution_sensitivity/ChIPseq_SET/sites/TFBS_separate.RData")
fit.SET_list = setfit
rm(setfit,petfit,exofit)

names(fit.exo_list) = names(peakList.exo)
names(fit.PET_list) = names(peakList.PET)
names(fit.SET_list) = names(peakList.SET)

seed = "12345"
rep = 2
maxdist = 400
mm = 3000

fit.exo = fit.exo_list[grep(seed,names(fit.exo_list))][[rep]]
fit.PET = fit.PET_list[grep(seed,names(fit.PET_list))][[1]]
fit.SET= fit.SET_list[grep(seed,names(fit.SET_list))][[1]]

peakList.exo = peakList.exo[grep(seed,names(peakList.exo))][[rep]]
peakList.PET = peakList.PET[grep(seed,names(peakList.PET))][[1]]
peakList.SET = peakList.SET[grep(seed,names(peakList.SET))][[1]]

# extract only coordinates from predictions

pred.exo <- data.table( "U00096", unlist(fit.exo@optMu) )
pred.PET <- data.table( "U00096", unlist(fit.PET@optMu) )
pred.SET <- data.table( "U00096", unlist(fit.SET@optMu) )

# IRanges

range.peak.exo = peakList.exo[,IRanges(start = V2,end = V3)]
range.peak.PET = peakList.PET[,IRanges(start = V2,end = V3)]
range.peak.SET = peakList.SET[,IRanges(start = V2,end = V3)]
range.prom = prom[,IRanges(start = V4,width = 1)]

# common peak region

overlap.prom.PET <- countOverlaps( range.peak.exo, range.peak.PET )
overlap.prom.SET <- countOverlaps( range.peak.exo, range.peak.SET )
table(overlap.prom.PET)
table(overlap.prom.SET)


## peakList.sub.exo <- peakList.exo[ overlap.prom.PET > 0 & overlap.prom.SET > 0 & V5 > 3000, ]

peakList.sub.exo <- peakList.exo[ overlap.prom.PET > 0 & overlap.prom.SET > 0 & V5 > mm, ]
nrow(peakList.sub.exo)

# match peak list & annotation

range.peak.exo <- peakList.sub.exo[,IRanges( start=V2, end=V3)]
overlap.prom.exo <-  countOverlaps( range.peak.exo, range.prom )
table( overlap.prom.exo )

peakList.sub <- peakList.sub.exo[ overlap.prom.exo > 0, ]
range.peak <- peakList.sub[,IRanges( start=V2, end=V3 )]

countVec <- countOverlaps( range.peak, range.prom )
matchMat <- as.matrix( findOverlaps( range.peak, range.prom ) )
matchInd <- split( matchMat[,2], matchMat[,1] )
table(countVec)

# resolution

result <- c()

for ( j in 1:length(countVec) ) {
    if (countVec[j] != 1 ) {
        next;
    }
    locp <- loc.prom[ matchInd[[ as.character(j) ]] ]
    
    result.j <- rep( NA, 3 )    
    if ( min(abs( locp - pred.exo[ !is.na(pred.exo$V2) ]$V2 )) <= 1000 ) {
        result.j[1] <- min(abs( locp - pred.exo[ !is.na(pred.exo$V2) ]$V2 ))
    }
    if ( min(abs( locp - pred.PET$V2 )) <= 1000 ) {
        result.j[2] <- min(abs( locp - pred.PET$V2 ))
    }
    if ( min(abs( locp - pred.SET$V2 )) <= 1000 ) {
        result.j[3] <- min(abs( locp - pred.SET$V2 ))
    }
    
    result <- rbind( result, result.j )
}

colnames(result) <- c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)" )

result_rene = result
row.names(result_rene) = NULL

result_rene = data.table(melt(result_rene))
result_rene[,Var1:=NULL]
setnames(result_rene,names(result_rene),c("Seq","Resolution"))

result_rene[, Seq := plyr::mapvalues(Seq,
   from = c("ChIP-exo","ChIP-Seq (PET)", "ChIP-Seq (SET)"),
   to = c("ChIP-exo","ChIP-Seq (PE)","ChIP-Seq (SE)"))]





## pdf( paste(dir.out,"resolution_exo_seq_R1.pdf",sep="") )
## #pdf( paste(dir.out,"resolution_exo_seq_R2.pdf",sep="") )
## boxplot( result, xlab="Dataset", ylab="Resolution", 
##     ylim=c(0,150), cex.lab=1.5, cex.main=1.5 )
## dev.off()

apply( result, 2, mean,trim = .1 )

# sensitivity

result <- c()
result.countVec <- c()

for ( j in 1:length(countVec) ) {
    if ( countVec[j] < 2 ) {
        next;
    }
    locp <- loc.prom[ matchInd[[ as.character(j) ]] ]
    #if ( min(abs(diff(locp))) < 50 ) {
    #    next;
    #}
    
    result.j <- matrix( NA, 3, length(locp) )  
    for ( k in 1:length(locp) ) {   
        result.j[ 1, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.exo[ !is.na(pred.exo$V2)]$V2 ) <= extension )) >= 1 )
        result.j[ 2, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.PET$V2 ) <= extension )) >= 1 )
        result.j[ 3, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.SET$V2 ) <= extension )) >= 1 )
    }
    
    result <- rbind( result, apply( result.j, 1, sum ) )
    result.countVec <- c( result.countVec, countVec[j] )
}

colnames(result) <- c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)" )
apply( result, 2, sum ) / sum(result.countVec)

# sensitivity

result <- c()
result.dist <- c()
result.min <- c()
result.countVec <- c()

for ( j in 1:length(countVec) ) {
    #if ( countVec[j] != 2 ) {
    #    next;
    #}
    #if ( countVec[j] < 1 ) {
	if ( countVec[j] < 2 ) {
        next;
    }
    locp <- loc.prom[ matchInd[[ as.character(j) ]] ]
    #if ( abs(diff(locp)) < 50 ) {
    #    next;
    #}
    
    result.j <- matrix( NA, 3, length(locp) )  
    dist.j <- matrix( NA, 3, length(locp) )
    for ( k in 1:length(locp) ) {   
        result.j[ 1, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.exo[ !is.na(pred.exo$V2) ]$V2 ) <= extension )) >= 1 )
        result.j[ 2, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.PET$V2 ) <= extension )) >= 1 )
        result.j[ 3, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.SET$V2) <= extension )) >= 1 )
        
        dist.j[ 1, k ] <- min(abs( locp[k] - pred.exo[ !is.na(pred.exo$V2)]$V2 ))
        dist.j[ 2, k ] <- min(abs( locp[k] - pred.PET$V2 ))
        dist.j[ 3, k ] <- min(abs( locp[k] - pred.SET$V2 ))
    }
    
    #result <- rbind( result, c( apply( result.j, 1, sum ), abs(diff(locp)) ) )
    result <- rbind( result, apply( result.j, 1, sum ) / length(locp) ) 
    result.dist <- c( result.dist, mean(abs(diff(sort(locp)))) )
    result.min <- rbind( result.min, apply( dist.j, 1, min ) )
    result.countVec <- c( result.countVec, countVec[j] )
}

colnames(result) <- c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)" )
apply( result, 2, mean )


result.final <- result[ 0 < result.dist & result.dist < maxdist, ]
result.dist.final <- result.dist[ 0 < result.dist & result.dist < maxdist ]
result.countVec.final <- result.countVec[ 0 < result.dist & result.dist < maxdist ]


dt =  data.table(melt(result.final))
dt[,Var1:= NULL]
setnames(dt,names(dt),c("Dataset","Sensitivity"))
dt[,distance:= result.dist.final]
dt[, Dataset := plyr::mapvalues(Dataset,
   from = c("ChIP-exo","ChIP-Seq (PET)", "ChIP-Seq (SET)"),
   to = c("exo","pet","set"))]


library(MASS)
library(dplyr)

dt  = dt %>%
  mutate(cate = ifelse(between(distance,0,100),
           "0 < distance <= 100",
           "distance > 100"))

dt = dt %>% mutate(cate = factor(cate,
      levels = c("0 < distance <= 100","distance > 100")))
  

dt  = dt %>%
  mutate(cate2 = ifelse(between(distance,0,100), "0 <= distance <= 100",
           ifelse(between(distance,100,200),"100 < distance <= 200",
                  "distance > 200")))

dt = dt %>% mutate(cate2 = factor(cate2,levels = c("0 <= distance <= 100",
                                        "100 < distance <= 200",
                                        "distance > 200")))


dt = data.table(dt)

dr = "figs/for_paper"
pdf(file = file.path(dr , paste0("subsample_resolution_old_data_",seed,"_rep2.pdf")),height = 3.5,width = 3.5)
p = ggplot(result_rene,aes(Seq,Resolution,colour = Seq))+
  geom_boxplot()+
  scale_y_continuous(limits = c(0,150))+xlab("")+
  theme_bw()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 6))+
  scale_color_brewer(palette = "Set1",name = "")+ggtitle("D")
u = print(p)
dev.off()





## ## pdf(file = file.path(dr,"sensitivity_exo_olda_data.pdf"),height = 3.5,width = 3.5)
dr = "figs/for_paper"
pdf(file = file.path(dr , paste0("subsample_sensitivity_old_data_",seed,"_rep2.pdf")),height = 3.5,width = 3.5)
p = ggplot(dt,aes(distance,Sensitivity,colour = Dataset,shape = Dataset))+
  geom_point(size = 2)+theme_bw()+
  stat_smooth(method = "nls",
              aes(group=Dataset),se =FALSE,size = 1,
              formula = y ~ a + b* x,
              method.args = list(
                algorithm = "port",
                start = c(a = 0.5,b = .1),
                lower = c(a = 0,b = 0),
                upper = c(a = 1,b = 1)) )+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0))+
  scale_color_brewer(palette = "Set1")+
  scale_shape_manual(values = c(1,3,4))+ggtitle("C")+
  xlab("Average distance between \nannotated binding events")
u = print(p)
dev.off()

## pdf(file = "Sensitivity_test2.pdf")
## p = ggplot(dt,aes(Dataset , Sensitivity,colour = Dataset))+
##   geom_boxplot()+geom_jitter()+
##   theme_bw()+
##   theme(legend.position = "top")+scale_color_brewer(palette = "Set1")+
##   facet_grid( ~ cate)
## u = print(p)
## dev.off()


## ggplot(dt,aes(Sensitivity,fill = Dataset))+
##   geom_histogram(aes(y=..density..),position = "dodge")+
##   theme_bw()+theme(legend.position = "top")+
##   scale_fill_brewer(palette = "Set1")+facet_grid( ~ cate)
  
  

## dev.off()

# COMPARISON with peak region (specificity)

## range.exo <- IRanges( start=pred.exo[ !is.na(pred.exo[,2]), 2 ], end=pred.exo[ !is.na(pred.exo[,2]), 2 ] )
## range.PET <- IRanges( start=pred.PET[,2], end=pred.PET[,2] )
## range.SET <- IRanges( start=pred.SET[,2], end=pred.SET[,2] )

## overlap.exo <- unlist( countOverlaps( range.peak, range.exo ) )
## overlap.PET <- unlist( countOverlaps( range.peak, range.PET ) )
## overlap.SET <- unlist( countOverlaps( range.peak, range.SET ) )

## table(overlap.exo)
## table(overlap.PET)
## table(overlap.SET)

## mean(overlap.exo)
## mean(overlap.PET)
## mean(overlap.SET)

## median(overlap.exo)
## median(overlap.PET)
## median(overlap.SET)

## overlap_mat <- cbind( overlap.exo, overlap.PET, overlap.SET )
## colnames(overlap_mat) <- c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)" )

## pdf( paste(dir.out,"npred_exo_seq_boxplot_R1.pdf",sep="") )
## boxplot( overlap_mat,
## 	xlab="Dataset", ylab="Average number of predictions",
## 	cex.lab=1.5, cex.axis=1.5, ylim=c(0,6) )
## dev.off()

## freq_overlap <- matrix( NA, 5, 3 )
## for ( i in 1:4 ) {
## 	freq_overlap[ i, 1 ] <- length(which( overlap.exo == (i-1) ))
## 	freq_overlap[ i, 2 ] <- length(which( overlap.PET == (i-1) ))
## 	freq_overlap[ i, 3 ] <- length(which( overlap.SET == (i-1) ))
## }
## freq_overlap[ 5, 1 ] <- length(which( overlap.exo >= 4 ))
## freq_overlap[ 5, 2 ] <- length(which( overlap.PET >= 4 ))
## freq_overlap[ 5, 3 ] <- length(which( overlap.SET >= 4 ))

## colnames(freq_overlap) <- c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)" )
## rownames(freq_overlap) <- c( "0", "1", "2", "3", ">= 4" )

## pdf( paste(dir.out,"npred_exo_seq_barplot_R1.pdf",sep="") )
## barplot( freq_overlap, beside=T, col=rainbow(5),
## 	xlab="Dataset", ylab="Frequency",
## 	cex.lab=1.5, cex.axis=1.5, ylim=c(0,200) )
## legend( "topleft", fill=rainbow(5), cex=1.3, bty="n",
## 	paste("# predictions",rownames(freq_overlap)) )
## dev.off()

