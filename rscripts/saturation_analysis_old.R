rm( list=ls() )
graphics.off()
library(mosaics)
library(IRanges)
library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# option

maxM <- 10
ext_iden <- 20
max_dist <- 1000
seedset <- c( 12345, 23456, 34567, 45678, 56789 )
ntop <- 500

dir_out <- "/p/keles/Dongjun/ChIP-exo/paper/figures/"

dir_mosaics_exo_R1 <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/mosaics_ChIPexo_R1/"
dir_mosaics_exo_R2 <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/mosaics_ChIPexo_R2/"
dir_mosaics_seq_PET <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/mosaics_ChIPseq_PET/"
dir_mosaics_seq_SET <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/mosaics_ChIPseq_SET/"

dir_dpeak_exo_R1 <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/dpeak_ChIPexo_R1/"
dir_dpeak_seq_PET <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/dpeak_ChIPseq_PET/"
dir_dpeak_seq_SET <- "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/dpeak_ChIPseq_SET/"

# read promoter regions for sigma70 from regulonDB

prom <- read.table(
    "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/database_regulonDB/PromoterSigma70Set.txt",
  skip=33, sep="\t")
loc.prom <- prom[,4]

# read peak summary

peaksummary_exo_R1 <- read.table( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/mosaics_ChIPexo_R1/peaksummary_subsample.txt", sep="\t" )
peaksummary_exo_R2 <- read.table( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/mosaics_ChIPexo_R2/peaksummary_subsample.txt", sep="\t" )
peaksummary_seq_PET <- read.table( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/mosaics_ChIPseq_PET/peaksummary_subsample.txt", sep="\t" )
peaksummary_seq_SET <- read.table( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/saturation_analysis/mosaics_ChIPseq_SET/peaksummary_subsample.txt", sep="\t" )

max_ss <- max( nrow(peaksummary_exo_R1), nrow(peaksummary_exo_R2), 
	nrow(peaksummary_seq_PET), nrow(peaksummary_seq_SET) )

# load peaks using all reads

load( paste(dir_mosaics_exo_R1,"peak_subsample_",nrow(peaksummary_exo_R1),".RData",sep="") )
peak_exo_R1 <- print(peak)
peak_exo_R1 <- peak_exo_R1[ order(peak_exo_R1$aveChipCount,decreasing=TRUE), ] #[ 1:1000, ]
rm( peak )

load( paste(dir_mosaics_exo_R2,"peak_subsample_",nrow(peaksummary_exo_R2),".RData",sep="") )
peak_exo_R2 <- print(peak)
peak_exo_R2 <- peak_exo_R2[ order(peak_exo_R2$aveChipCount,decreasing=TRUE), ] #[ 1:1000, ]
rm( peak )

load( paste(dir_mosaics_seq_PET,"peak_subsample_",nrow(peaksummary_seq_PET),".RData",sep="") )
peak_seq_PET <- print(peak)
peak_seq_PET <- peak_seq_PET[ order(peak_seq_PET$aveChipCount,decreasing=TRUE), ] #[ 1:1000, ]
rm( peak )

load( paste(dir_mosaics_seq_SET,"peak_subsample_",nrow(peaksummary_seq_SET),".RData",sep="") )
peak_seq_SET <- print(peak)
peak_seq_SET <- peak_seq_SET[ order(peak_seq_SET$aveChipCount,decreasing=TRUE), ] #[ 1:1000, ]
rm( peak )

# extract common peaks

range_exo_R1 <- IRanges( start=peak_exo_R1[,2], end=peak_exo_R1[,3] )
range_exo_R2 <- IRanges( start=peak_exo_R2[,2], end=peak_exo_R2[,3] )
range_seq_PET <- IRanges( start=peak_seq_PET[,2], end=peak_seq_PET[,3] )
range_seq_SET <- IRanges( start=peak_seq_SET[,2], end=peak_seq_SET[,3] )

loc_common <- which( countOverlaps( range_exo_R1, range_exo_R2 ) > 0 &
	countOverlaps( range_exo_R1, range_seq_PET ) > 0 &
	countOverlaps( range_exo_R1, range_seq_SET ) > 0 )

peak_common <- peak_exo_R1[ loc_common, ]

# ground truth: events

range_common <- IRanges( start=peak_common[,2], end=peak_common[,3] )
range_prom <- IRanges( start=loc.prom, end=loc.prom )
truth_prom <- loc.prom[ countOverlaps( range_prom, range_common ) > 0 ]
	
# overlap with common peaks

npred <- npeak <- iden <- resol <- matrix( NA, maxM, 3 )
iden.all <- rep( NA, maxM )
range_target <- IRanges( start=truth_prom - ext_iden, end=truth_prom + ext_iden )

for ( i in 1:10 ) {
	print(i)
	
	if ( i < 10 ) {
		pid <- i + 1
	} else {
		pid <- 1
	}
	
	# load mosaics peak lists
	
	load( paste(dir_mosaics_exo_R1,"peak_0.1M_seed_12345_subsample_1.RData",sep="") )
	peak_exo_R1 <- print(peak)
	rm( peak )
	
	load( paste(dir_mosaics_seq_PET,"peak_0.1M_seed_12345_subsample_1.RData",sep="") )
	peak_seq_PET <- print(peak)
	rm( peak )
	
	load( paste(dir_mosaics_seq_SET,"peak_0.1M_seed_12345_subsample_1.RData",sep="") )
	peak_seq_SET <- print(peak)
	rm( peak )
	
	# determine top peaks
	
	peak_exo_R1 <- peak_exo_R1[ order( peak_exo_R1$aveChipCount, decreasing=TRUE ), ]
	peak_exo_R1 <- peak_exo_R1[ 1:500, ]
	
	peak_seq_PET <- peak_seq_PET[ order( peak_seq_PET$aveChipCount, decreasing=TRUE ), ]
	peak_seq_PET <- peak_seq_PET[ 1:500, ]
	
	peak_seq_SET <- peak_seq_SET[ order( peak_exo_R1$aveChipCount, decreasing=TRUE ), ]
	peak_seq_SET <- peak_seq_SET[ 1:500, ]
	
	crit_exo_R1 <- sort( peak_exo_R1$aveChipCount, decreasing=TRUE )[ntop]
	ind_exo_R1 <- as.numeric( peak_exo_R1$aveChipCount >= crit_exo_R1 )
	
	crit_seq_PET <- sort( peak_seq_PET$aveChipCount, decreasing=TRUE )[ntop]
	ind_seq_PET <- as.numeric( peak_seq_PET$aveChipCount >= crit_seq_PET )
	
	crit_seq_SET <- sort( peak_seq_SET$aveChipCount, decreasing=TRUE )[ntop]
	ind_seq_SET <- as.numeric( peak_seq_SET$aveChipCount >= crit_seq_SET )
	
	# summary over multiple seeds
	
	npeak_mat <- npred_mat <- iden_mat <- resol_mat <- matrix( NA, length(seedset), 3 )
	iden_all_mat <- rep( NA, length(seedset) )
	
	for ( ss in 1:length(seedset) ) {
	
		seedj <- seedset[ss]
		
		# load dpeak predictions
		
		load( file=paste(dir_dpeak_exo_R1,"fit_0.01M_seed_",seedj,"_subsample_",pid,".RData",sep="") )
		optMu <- fitSET@optMu
		npeak.i <- 0
		pred_exo_R1 <- c()
		for ( j in 1:length(optMu) ) {
			if ( ind_exo_R1[j] ==1 & !is.na(optMu[[j]][1]) ) {
				npeak.i <- npeak.i + 1
				pred_exo_R1 <- c( pred_exo_R1, optMu[[j]] )
			}
		}
		#npeak[ i, 1 ] <- npeak.i	
		npeak_mat[ ss, 1 ] <- npeak.i
		pred_exo_R1 <- sort(pred_exo_R1)
		rm(fitSET)
		
		load( file=paste(dir_dpeak_seq_PET,"fit_0.01M_seed_",seedj,"_subsample_",pid,".RData",sep="") )
		optMu <- fitPET@optMu
		npeak.i <- 0
		pred_seq_PET <- c()
		for ( j in 1:length(optMu) ) {
			if ( ind_seq_PET[j] ==1 & !is.na(optMu[[j]][1]) ) {
				npeak.i <- npeak.i + 1
				pred_seq_PET <- c( pred_seq_PET, optMu[[j]] )
			}
		}
		#npeak[ i, 2 ] <- npeak.i	
		npeak_mat[ ss, 2 ] <- npeak.i
		pred_seq_PET <- sort(pred_seq_PET)
		rm(fitPET)
		
		load( file=paste(dir_dpeak_seq_SET,"fit_0.01M_seed_",seedj,"_subsample_",pid,".RData",sep="") )
		optMu <- fitSET@optMu
		npeak.i <- 0
		pred_seq_SET <- c()
		for ( j in 1:length(optMu) ) {
			if ( ind_seq_SET[j] ==1 & !is.na(optMu[[j]][1]) ) {
				npeak.i <- npeak.i + 1
				pred_seq_SET <- c( pred_seq_SET, optMu[[j]] )
			}
		}
		#npeak[ i, 3 ] <- npeak.i	
		npeak_mat[ ss, 3 ] <- npeak.i
		pred_seq_SET <- sort(pred_seq_SET)
		rm(fitSET)
		
		npred_mat[ ss, 1 ] <- length(pred_exo_R1)
		npred_mat[ ss, 2 ] <- length(pred_seq_PET)
		npred_mat[ ss, 3 ] <- length(pred_seq_SET)
		
		# number of identified targets
		
		range_exo_R1 <- IRanges( start=pred_exo_R1, end=pred_exo_R1 )
		range_seq_PET <- IRanges( start=pred_seq_PET, end=pred_seq_PET )
		range_seq_SET <- IRanges( start=pred_seq_SET, end=pred_seq_SET )
		
		iden_mat[ ss, 1 ] <- sum( countOverlaps( range_target, range_exo_R1 ) > 0 )
		iden_mat[ ss, 2 ] <- sum( countOverlaps( range_target, range_seq_PET ) > 0 )
		iden_mat[ ss, 3 ] <- sum( countOverlaps( range_target, range_seq_SET ) > 0 )
		
		pred_all <- c( pred_exo_R1, pred_seq_PET, pred_seq_SET )
		#pred_all <- c( pred_exo_R1, pred_seq_PET )
		range_all <- IRanges( start=pred_all, end=pred_all )
		iden_all_mat[ ss ] <- sum( countOverlaps( range_target, range_all ) > 0 )
		
		# resolution
		
		dist_exo_R1 <- dist_seq_PET <- dist_seq_SET <- rep( NA, length(truth_prom) )
		
		for ( j in 1:length(truth_prom) ) {
			if ( min(abs( pred_exo_R1 - truth_prom[j] )) <= max_dist ) {
				dist_exo_R1[j] <- min(abs( pred_exo_R1 - truth_prom[j] ))
			}
			if ( min(abs( pred_seq_PET - truth_prom[j] )) <= max_dist ) {
				dist_seq_PET[j] <- min(abs( pred_seq_PET - truth_prom[j] ))
			}
			if ( min(abs( pred_seq_SET - truth_prom[j] )) <= max_dist ) {
				dist_seq_SET[j] <- min(abs( pred_seq_SET - truth_prom[j] ))
			}
		}
		
		resol_mat[ ss, 1 ] <- median( dist_exo_R1, na.rm=TRUE )
		resol_mat[ ss, 2 ] <- median( dist_seq_PET, na.rm=TRUE )
		resol_mat[ ss, 3 ] <- median( dist_seq_SET, na.rm=TRUE )
	}
	
	# update summary statistics
		
	npeak[ i,  ] <- apply( npeak_mat, 2, mean )
	npred[ i,  ] <- apply( npred_mat, 2, mean )
	iden[ i,  ] <- apply( iden_mat, 2, mean )
	resol[ i,  ] <- apply( resol_mat, 2, mean )
	iden.all[ i ] <- mean(iden_all_mat)
}

# plot number of candidate regions

dr = "/p/keles/ChIPexo/volume3/ChIPexo/figs/for_paper"

what <- c( "exo", "pet", "set")

dt1 = data.table(npeak)
setnames(dt1,names(dt1),what)
dt1 = melt(dt1)
dt1[,nreads:= 0.01 * 1:maxM]
setnames(dt1,names(dt1),c("Dataset","candre","nreads"))

dt2 = data.table(npred)
setnames(dt2,names(dt2),what)
dt2 = melt(dt2)
dt2[,nreads:= 0.01 * 1:maxM]
setnames(dt2,names(dt2),c("Dataset","predre","nreads"))

dt3 = data.table(iden)# / max(iden.all))
setnames(dt3,names(dt3),what)
dt3 = melt(dt3)
dt3[,nreads:= 0.01 * 1:maxM]
setnames(dt3,names(dt3),c("Dataset","iden","nreads"))

dt4 = data.table(resol)
setnames(dt4,names(dt4),what)
dt4 = melt(dt4)
dt4[,nreads:= 0.01 * 1:maxM]
setnames(dt4,names(dt4),c("Dataset","resol","nreads"))

library(grid)
library(gridExtra)
library(RColorBrewer)

pdf(file = "figs/for_paper/saturation_analysis_old.pdf")
p1 <- ggplot(dt1,aes(nreads,candre, colour = Dataset))+
  geom_line(size = 1.2)+theme_bw()+scale_color_brewer(palette = "Set1",name="")+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Number of reads (million)")+ylab("Number of candidate regions")+
  ylim(0,600)+ggtitle("A")
p2 <- ggplot(dt2,aes(nreads,predre, colour = Dataset))+
  geom_line(size = 1.2)+theme_bw()+scale_color_brewer(palette = "Set1",name="")+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Number of reads (million)")+ylab("Number of predicted events")+
  ylim(0,1500)+ggtitle("B")
p3 <- ggplot(dt3,aes(nreads,iden, colour = Dataset))+
  geom_line(size = 1.2)+theme_bw()+scale_color_brewer(palette = "Set1",name="")+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Number of reads (million)")+ylab("Number of identified targets")+
  ylim(0,250)+ggtitle("C")
p4 <- ggplot(dt4,aes(nreads,resol, colour = Dataset))+
  geom_line(size = 1.2)+theme_bw()+scale_color_brewer(palette = "Set1",name="")+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Number of reads (million)")+ylab("Resolution")+
  ylim(0,50)+ggtitle("D")
grid.arrange(p1,p2,p3,p4,nrow =2)
dev.off()



pdf( paste(dir_out,"saturation_regulonDB_ncand_ntop_",ntop,".pdf",sep="") )

plot( 0.01 * c(1:maxM), 0.01 * c(1:maxM), xlim=0.01 * c(1,maxM), ylim=c(0,600), type="n",
	xlab="Number of reads (million)", ylab="Number of candidate regions",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5 )
lines( 0.01 * c(1:maxM), npeak[,1], col=1, lty=1, lwd=3 )
lines( 0.01 * c(1:maxM), npeak[,2], col=2, lty=2, lwd=3 )
lines( 0.01 * c(1:maxM), npeak[,3], col="blue", lty=4, lwd=3 )
legend( "bottomright", col=c(1,2,"blue"), lty=c(1,2,4), lwd=3, cex=1.5,
	c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)") )

dev.off()

# plot number of predictions

pdf( paste(dir_out,"saturation_regulonDB_npred_ntop_",ntop,".pdf",sep="") )
plot( 0.01 * c(1:maxM), 0.01 * c(1:maxM), xlim=0.01 * c(1,maxM), ylim=c(0,1500), type="n",
	xlab="Number of reads (million)", ylab="Number of predicted events",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5 )
lines( 0.01 * c(1:maxM), npred[,1], col=1, lty=1, lwd=3 )
lines( 0.01 * c(1:maxM), npred[,2], col=2, lty=2, lwd=3 )
lines( 0.01 * c(1:maxM), npred[,3], col="blue", lty=4, lwd=3 )
legend( "bottomright", col=c(1,2,"blue"), lty=c(1,2,4), lwd=3, cex=1.5,
	c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)") )
dev.off()
	
# plot number of identified targets

pdf( paste(dir_out,"saturation_regulonDB_propiden_ntop_",ntop,".pdf",sep="") )
plot( 0.01 * c(1:maxM), 0.01 * c(1:maxM), xlim=0.01 * c(1,maxM), ylim=c(0,1), type="n",
	xlab="Number of reads (million)", ylab="Proportion of identified targets",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5 )
lines( 0.01 * c(1:maxM), iden[,1] / max(iden.all), col=1, lty=1, lwd=3 )
lines( 0.01 * c(1:maxM), iden[,2] / max(iden.all), col=2, lty=2, lwd=3 )
lines( 0.01 * c(1:maxM), iden[,3] / max(iden.all), col="blue", lty=4, lwd=3 )
#abline( h=length(truth_prom) )
legend( "bottomright", col=c(1,2,"blue"), lty=c(1,2,4), lwd=3, cex=1.5,
	c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)") )
dev.off()

# plot resolution

pdf( paste(dir_out,"saturation_regulonDB_resol_ntop_",ntop,".pdf",sep="") )
plot( 0.01 * c(1:maxM), 0.01 * c(1:maxM), xlim=0.01 * c(1,maxM), ylim=c(0,50), type="n",
	xlab="Number of reads (million)", ylab="Resolution",
	cex.lab=1.5, cex.axis=1.5, cex.main=1.5 )
lines( 0.01 * c(1:maxM), resol[,1], col=1, lty=1, lwd=3 )
lines( 0.01 * c(1:maxM), resol[,2], col=2, lty=2, lwd=3 )
lines( 0.01 * c(1:maxM), resol[,3], col="blue", lty=4, lwd=3 )
legend( "topright", col=c(1,2,"blue"), lty=c(1,2,4), lwd=3, cex=1.5,
	c( "ChIP-exo", "ChIP-Seq (PET)", "ChIP-Seq (SET)") )
dev.off()
