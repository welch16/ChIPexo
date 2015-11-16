rm( list=ls() )
graphics.off()
library(dpeak)
library(IRanges)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# options

extension <- 20

#dir.peak <- "/scratch/deconvolution/comparison_qSET/"
dir.out <- "/p/keles/Dongjun/ChIP-exo/paper/figures/"

# read promoter regions for sigma70 from regulonDB

prom <- read.table(
    "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/database_regulonDB/PromoterSigma70Set.txt",
  skip=33, sep="\t")
loc.prom <- prom[,4]

# load peak region & motif

peakList.exo <- read.table( 
    "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_mosaics/ChIP-exo_sigma70_exp_phase_R1_mosaics_1S_peak.bed",
    sep="\t", stringsAsFactors=FALSE )

peakList.PET <- read.table( 
    "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/chipseq/Kiley_PET/sigma70+O2/mosaics_PET/peaklist_IO.bed",
    sep="\t", stringsAsFactors=FALSE, skip=1 )

peakList.SET <- read.table( 
    "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/chipseq/Kiley_SET/sigma70+O2/mosaics_SET/peaklist_IO.bed",
    sep="\t", stringsAsFactors=FALSE, skip=1 )

# load predictions

load( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_R1/ChIP-exo_sigma70_exp_phase_dpeak_fit.RData" )
#load( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_R2/ChIP-exo_sigma70_exp_phase_dpeak_fit.RData" )
fit.exo <- fitSET
rm(fitSET)

load( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/chipseq/Kiley_PET/sigma70+O2/mosaics_PET/sigma70+O2_PET_fit.RData" )
fit.PET <- fitPET
rm(fitPET)

load( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/chipseq/Kiley_SET/sigma70+O2/mosaics_SET/fit_sigma70+O2.RData" )
fit.SET <- fitSET
rm(fitSET)
    
# extract only coordinates from predictions

pred.exo <- data.frame( "U00096", unlist(fit.exo@optMu) )
pred.PET <- data.frame( "U00096", unlist(fit.PET@optMu) )
pred.SET <- data.frame( "U00096", unlist(fit.SET@optMu) )

# IRanges

range.peak.exo <- IRanges( start=peakList.exo[,2], end=peakList.exo[,3] )
range.peak.PET <- IRanges( start=peakList.PET[,2], end=peakList.PET[,3] )
range.peak.SET <- IRanges( start=peakList.SET[,2], end=peakList.SET[,3] )
range.prom <- IRanges( start=prom[,4], end=prom[,4] )

# common peak region

overlap.prom.PET <- countOverlaps( range.peak.exo, range.peak.PET )
overlap.prom.SET <- countOverlaps( range.peak.exo, range.peak.SET )
table(overlap.prom.PET)
table(overlap.prom.SET)

#peakList.sub.exo <- peakList.exo[ overlap.prom.PET > 0 & overlap.prom.SET > 0, ]
peakList.sub.exo <- peakList.exo[ overlap.prom.PET > 0 & overlap.prom.SET > 0 &
	peakList.exo[,5] > 3000, ]
nrow(peakList.sub.exo)

# match peak list & annotation

range.peak.exo <- IRanges( start=peakList.sub.exo[,2], end=peakList.sub.exo[,3] )
overlap.prom.exo <-  countOverlaps( range.peak.exo, range.prom )
table( overlap.prom.exo )

peakList.sub <- peakList.sub.exo[ overlap.prom.exo > 0, ]
range.peak <- IRanges( start=peakList.sub[,2], end=peakList.sub[,3] )

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
    if ( min(abs( locp - pred.exo[ !is.na(pred.exo[,2]), 2 ] )) <= 1000 ) {
        result.j[1] <- min(abs( locp - pred.exo[ !is.na(pred.exo[,2]), 2 ] ))
    }
    if ( min(abs( locp - pred.PET[,2] )) <= 1000 ) {
        result.j[2] <- min(abs( locp - pred.PET[,2] ))
    }
    if ( min(abs( locp - pred.SET[,2] )) <= 1000 ) {
        result.j[3] <- min(abs( locp - pred.SET[,2] ))
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
   to = c("exo","pet","set"))]


dr = "figs/for_paper"
pdf(file = file.path(dr , "resolution_by_dataset_old_data.pdf"),height = 3.5,width = 3.5)
p = ggplot(result_rene,aes(Seq,Resolution,colour = Seq))+geom_boxplot()+
  scale_y_continuous(limits = c(0,150))+xlab("")+
  theme_bw()+
  theme(legend.position = "none",plot.title = element_text(hjust = 0))+
  scale_color_brewer(palette = "Set1")+ggtitle("B")
u = print(p)
dev.off()


pdf( paste(dir.out,"resolution_exo_seq_R1.pdf",sep="") )
#pdf( paste(dir.out,"resolution_exo_seq_R2.pdf",sep="") )
boxplot( result, xlab="Dataset", ylab="Resolution", 
    ylim=c(0,150), cex.lab=1.5, cex.main=1.5 )
dev.off()

apply( result, 2, mean )

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
            as.numeric( length(which( abs( locp[k] - pred.exo[ !is.na(pred.exo[,2]) ,2 ] ) <= extension )) >= 1 )
        result.j[ 2, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.PET[,2] ) <= extension )) >= 1 )
        result.j[ 3, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.SET[,2] ) <= extension )) >= 1 )
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
            as.numeric( length(which( abs( locp[k] - pred.exo[ !is.na(pred.exo[,2]) ,2 ] ) <= extension )) >= 1 )
        result.j[ 2, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.PET[,2] ) <= extension )) >= 1 )
        result.j[ 3, k ] <- 
            as.numeric( length(which( abs( locp[k] - pred.SET[,2] ) <= extension )) >= 1 )
        
        dist.j[ 1, k ] <- min(abs( locp[k] - pred.exo[ !is.na(pred.exo[,2]), 2 ] ))
        dist.j[ 2, k ] <- min(abs( locp[k] - pred.PET[,2] ))
        dist.j[ 3, k ] <- min(abs( locp[k] - pred.SET[,2] ))
    }
    
    #result <- rbind( result, c( apply( result.j, 1, sum ), abs(diff(locp)) ) )
    result <- rbind( result, apply( result.j, 1, sum ) / length(locp) ) 
    result.dist <- c( result.dist, mean(abs(diff(sort(locp)))) )
    result.min <- rbind( result.min, apply( dist.j, 1, min ) )
    result.countVec <- c( result.countVec, countVec[j] )
}

colnames(result) <- c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)" )
apply( result, 2, mean )

maxdist <- 450

result.final <- result[ 0 < result.dist & result.dist < maxdist, ]
result.dist.final <- result.dist[ 0 < result.dist & result.dist < maxdist ]
result.countVec.final <- result.countVec[ 0 < result.dist & result.dist < maxdist ]

pdf( paste(dir.out,"sensitivity_exo_seq_R1.pdf",sep="") )
#pdf( paste(dir.out,"sensitivity_exo_seq_R2.pdf",sep="") )

library(MASS)
plot( 1:5, 1:5, 
	xlim=c(min(result.dist.final),max(result.dist.final)), ylim=c(0,1), type="n",
	xlab="Average distance between annotated binding events", ylab="Sensitivity",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5 )
points( result.dist.final, result.final[,1], col="gray", pch=1, cex=2 )
points( result.dist.final, result.final[,2], col="pink", pch=3, cex=2 )
points( result.dist.final, result.final[,3], col="lightblue", pch=4, cex=2 )
abline( rlm( result.final[,1] ~ result.dist.final ), col="black", lwd=3 )
abline( rlm( result.final[,2] ~ result.dist.final ), col="red", lwd=3 )
abline( rlm( result.final[,3] ~ result.dist.final ), col="blue", lwd=3 )
legend( "topright", lwd=3, col=c("black","red","blue"), 
	colnames(result), cex=1.5, bg="white" )

dev.off()

dt =  data.table(melt(result.final))
dt[,Var1:= NULL]
setnames(dt,names(dt),c("Dataset","Sensitivity"))
dt[,distance:= result.dist.final]
dt[, Dataset := plyr::mapvalues(Dataset,
   from = c("ChIP-exo","ChIP-Seq (PET)", "ChIP-Seq (SET)"),
   to = c("exo","pet","set"))]


library(MASS)

pdf(file = file.path(dr,"sensitivity_exo_olda_data.pdf"),height = 3.5,width = 3.5)
p = ggplot(dt,aes(distance,Sensitivity,colour = Dataset,shape = Dataset))+
  geom_point(alpha = I(1/2),size = 2)+theme_bw()+
  geom_smooth(method = "rlm",aes(group=Dataset),se =FALSE,size = 1)+
  theme(legend.position = "none",plot.title = element_text(hjust = 0))+scale_color_brewer(palette = "Set1")+
  scale_shape_manual(values = c(1,3,4))+ggtitle("A")+
  xlab("Average distance between \nannotated binding events")
u = print(p)
dev.off()

# COMPARISON with peak region (specificity)

range.exo <- IRanges( start=pred.exo[ !is.na(pred.exo[,2]), 2 ], end=pred.exo[ !is.na(pred.exo[,2]), 2 ] )
range.PET <- IRanges( start=pred.PET[,2], end=pred.PET[,2] )
range.SET <- IRanges( start=pred.SET[,2], end=pred.SET[,2] )

overlap.exo <- unlist( countOverlaps( range.peak, range.exo ) )
overlap.PET <- unlist( countOverlaps( range.peak, range.PET ) )
overlap.SET <- unlist( countOverlaps( range.peak, range.SET ) )

table(overlap.exo)
table(overlap.PET)
table(overlap.SET)

mean(overlap.exo)
mean(overlap.PET)
mean(overlap.SET)

median(overlap.exo)
median(overlap.PET)
median(overlap.SET)

overlap_mat <- cbind( overlap.exo, overlap.PET, overlap.SET )
colnames(overlap_mat) <- c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)" )

pdf( paste(dir.out,"npred_exo_seq_boxplot_R1.pdf",sep="") )
boxplot( overlap_mat,
	xlab="Dataset", ylab="Average number of predictions",
	cex.lab=1.5, cex.axis=1.5, ylim=c(0,6) )
dev.off()

freq_overlap <- matrix( NA, 5, 3 )
for ( i in 1:4 ) {
	freq_overlap[ i, 1 ] <- length(which( overlap.exo == (i-1) ))
	freq_overlap[ i, 2 ] <- length(which( overlap.PET == (i-1) ))
	freq_overlap[ i, 3 ] <- length(which( overlap.SET == (i-1) ))
}
freq_overlap[ 5, 1 ] <- length(which( overlap.exo >= 4 ))
freq_overlap[ 5, 2 ] <- length(which( overlap.PET >= 4 ))
freq_overlap[ 5, 3 ] <- length(which( overlap.SET >= 4 ))

colnames(freq_overlap) <- c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)" )
rownames(freq_overlap) <- c( "0", "1", "2", "3", ">= 4" )

pdf( paste(dir.out,"npred_exo_seq_barplot_R1.pdf",sep="") )
barplot( freq_overlap, beside=T, col=rainbow(5),
	xlab="Dataset", ylab="Frequency",
	cex.lab=1.5, cex.axis=1.5, ylim=c(0,200) )
legend( "topleft", fill=rainbow(5), cex=1.3, bty="n",
	paste("# predictions",rownames(freq_overlap)) )
dev.off()

