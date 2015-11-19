rm( list=ls() )
graphics.off()
library(mosaics)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(reshape2)

Rlist <- system( "ls /u/w/e/welch/Desktop/Docs/packages/mosaics/R/*.R", intern=TRUE )
for ( i in 1:length(Rlist) ) {
	source( Rlist[i] )
}


#constructBins( infile="CTCF.bowtie", fileFormat="bowtie", outfileLoc="./", 
#	byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
#	PET=FALSE, fragLen=200, binSize=200, capping=0, perl = "perl" )
         
dir_out <- "/p/keles/Dongjun/ChIP-exo/paper/figures/"

bin.exo <- readBins( c("chip","M","GC","N"),
	c( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/ChIP-exo/CTCF.bowtie_fragL200_bin200.txt",
	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/map50/map50_fragL200_bin200.txt",
	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/GC/GC_fragL200_bin200.txt",
	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/N/N_fragL200_bin200.txt" ) )

bin.seq <- readBins( c("chip","input","M","GC","N"),
	c( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/ChIP-seq_Crawford/CTCF_Crawford.bowtie_fragL200_bin200.txt",
	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/Input_Crawford/Input_Crawford.bowtie_fragL200_bin200.txt",
	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/map50/map50_fragL200_bin200.txt",
	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/GC/GC_fragL200_bin200.txt",
	"/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/CTCF/raw_data/N/N_fragL200_bin200.txt" ) )

statM <- .computeStat( Y=tagCount(bin.exo), S=mappability(bin.exo) ) 
statGC <- .computeStat( Y=tagCount(bin.exo), S=gcContent(bin.exo) )

M <- data.table( mapp = statM$uS, mean = statM$meanYall,
  Var = statM$varYall , N = statM$nitem)
GC <- data.table( gc = statGC$uS, mean = statGC$meanYall,
  Var = statGC$varYall , N = statGC$nitem)

M[,lb := mean - 1.96 * sqrt(Var / N)]
M[,ub := mean + 1.96 * sqrt(Var / N)]

GC[,lb := mean - 1.96 * sqrt(Var / N)]
GC[,ub := mean + 1.96 * sqrt(Var / N)]


  ## lb = statM$meanYall-1.96*sqrt(statM$varYall/statM$nitem),
  ## ub = statM$meanYall+1.96*sqrt(statM$varYall/statM$nitem))

pdf(file = "figs/for_paper/eukaryotic_bias_CTCF.pdf",width = 4 , height = 4)
ggplot( M , aes( x = mapp ,ymin = lb , ymax = ub))+geom_linerange()+xlim(0,1)+
  geom_point(aes( x = mapp, y = mean),size = 2)+theme_bw()+
  theme(plot.title = element_text(hjust = 0))+ggtitle("C")+
  coord_cartesian(xlim = c(0,1), ylim = c(-.1,7))+
  xlab("Mappability score")+ylab("Mean ChIP tag count")
ggplot( GC , aes( x = gc ,ymin = lb , ymax = ub))+geom_linerange()+xlim(0,1)+
  geom_point(aes( x =gc, y = mean),size = 2)+theme_bw()+
  theme(plot.title = element_text(hjust = 0))+ggtitle("D")+
  coord_cartesian(xlim = c(0,1), ylim = c(-.2,15))+
  xlab("GC content score")+ylab("Mean ChIP tag count")
dev.off()




## pdf( paste(dir_out,"seq_bias_CTCF_mappability.pdf",sep="") )        
## plot( statM$uS, statM$meanYall,
##     xlab='Mappability score', ylab='Mean ChIP tag count', 
##     #main='Mappability score vs. Mean ChIP tag count',
##     ylim=quantile( statM$meanYall, prob=c(0.05,0.95) ),
##     cex.axis=1.5, cex.lab=1.5 )
## segments( statM$uS, statM$meanYall, 
##     statM$uS, statM$meanYall+1.96*sqrt(statM$varYall/statM$nitem) )
## segments( statM$uS, statM$meanYall, 
##     statM$uS, statM$meanYall-1.96*sqrt(statM$varYall/statM$nitem) )
## dev.off()

## pdf( paste(dir_out,"seq_bias_CTCF_GC.pdf",sep="") )
## plot( statGC$uS, statGC$meanYall,
##     xlab='GC content score', ylab='Mean ChIP tag count', 
##     #main='GC content score vs. Mean ChIP tag count',
##     ylim=quantile( statGC$meanYall, prob=c(0.05,0.95) ),
##     cex.axis=1.5, cex.lab=1.5 )
## segments( statGC$uS, statGC$meanYall, 
##     statGC$uS, statGC$meanYall+1.96*sqrt(statGC$varYall/statGC$nitem) )
## segments( statGC$uS, statGC$meanYall, 
##     statGC$uS, statGC$meanYall-1.96*sqrt(statGC$varYall/statGC$nitem) )
## dev.off()
## #plot( bin.seq, plotType="M" )
## #plot( bin.seq, plotType="GC" )
	
