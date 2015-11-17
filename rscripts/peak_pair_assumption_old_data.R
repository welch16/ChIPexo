rm( list=ls() )
graphics.off()
library(dpeak)
library(hexbin)
library(IRanges)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

dir.out <- "/p/keles/Dongjun/ChIP-exo/paper/figures/"

# load peak lists

peak.R1 <- read.table( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_mosaics/ChIP-exo_sigma70_exp_phase_R1_mosaics_1S_peak.bed", sep="\t" )
peak.R2 <- read.table( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_mosaics/ChIP-exo_sigma70_exp_phase_R2_mosaics_1S_peak.bed", sep="\t" )

# load predictions

load( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_R1/ChIP-exo_sigma70_exp_phase_dpeak_fit.RData" )
fit.exo1 <- fitSET
rm(fitSET)
load( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/analysis_sigma70/analysis_sigma70_R2/ChIP-exo_sigma70_exp_phase_dpeak_fit.RData" )
fit.exo2 <- fitSET
rm(fitSET)

load( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/chipseq/Kiley_PET/sigma70+O2/mosaics_PET/sigma70+O2_PET_fit.RData" )
fit.PET <- fitPET
rm(fitPET)

load( "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/chipseq/Kiley_SET/sigma70+O2/mosaics_SET/fit_sigma70+O2.RData" )
fit.SET <- fitSET
rm(fitSET)

# extract delta & sigma

delta.exo1 <- unlist(fit.exo1@optDelta)
sigma.exo1 <- unlist(fit.exo1@optSigma)

delta.exo2 <- unlist(fit.exo2@optDelta)
sigma.exo2 <- unlist(fit.exo2@optSigma)

delta.SET <- unlist(fit.SET@optDelta)
sigma.SET <- unlist(fit.SET@optSigma)

# strength

str1 <- sapply( fit.exo1@fragSet, nrow )
str1 <- str1[ !is.na(delta.exo1) ]

str2 <- sapply( fit.exo2@fragSet, nrow )
str2 <- str2[ !is.na(delta.exo2) ]

# strand balance

balance1 <- sapply( fit.exo1@fragSet, function(x) 
	if ( !is.na(x[1,1]) ) {
		return( mean( x[,3]=="F" ) )
	} else {
		return( NA )
	} )
balance1 <- balance1[ !is.na(delta.exo1) ]

balance2 <- sapply( fit.exo2@fragSet, function(x) 
	if ( !is.na(x[1,1]) ) {
		return( mean( x[,3]=="F" ) )
	} else {
		return( NA )
	} )
balance2 <- balance2[ !is.na(delta.exo2) ]

# delta & sigma

delta.exo1 <- delta.exo1[ !is.na(delta.exo1) ]
sigma.exo1 <- sigma.exo1[ !is.na(sigma.exo1) ]

delta.exo2 <- delta.exo2[ !is.na(delta.exo2) ]
sigma.exo2 <- sigma.exo2[ !is.na(sigma.exo2) ]

delta.SET <- delta.SET[ !is.na(delta.SET) ]
sigma.SET <- sigma.SET[ !is.na(sigma.SET) ]

# plot

pdf( paste(dir.out,"issue_delta_estimate_exo_R1.pdf",sep="") )

a <- hexbin( delta.exo1, sigma.exo1, xbins = 30)
hbp <- plot( a, trans = log, inv = exp, 
	xlab = "delta estimate", ylab = "sigma estimate", 
	colramp = rainbow)

plot( delta.exo1, str1, log="y",
	xlab="delta estimate", ylab="Number of reads in the region",
	cex.axis=1.5, cex.lab=1.5 )

plot( sigma.exo1, str1, log="y",
	xlab="sigma estimate", ylab="Number of reads in the region",
	cex.axis=1.5, cex.lab=1.5 )

plot( delta.exo1, balance1, 
	xlab="delta estimate", ylab="Proportion of forward strand reads",
	cex.axis=1.5, cex.lab=1.5 )

plot( sigma.exo1, balance1, 
	xlab="sigma estimate", ylab="Proportion of forward strand reads",
	cex.axis=1.5, cex.lab=1.5 )

dev.off()

pdf( paste(dir.out,"issue_delta_estimate_exo_R2.pdf",sep="") )

a <- hexbin( delta.exo2, sigma.exo2, xbins = 30)
hbp <- plot( a, trans = log, inv = exp, 
	xlab = "delta estimate", ylab = "sigma estimate", 
	colramp = rainbow)

plot( delta.exo2, str2, log="y",
	xlab="delta estimate", ylab="Number of reads in the region",
	cex.axis=1.5, cex.lab=1.5 )

plot( sigma.exo2, str2, log="y",
	xlab="sigma estimate", ylab="Number of reads in the region",
	cex.axis=1.5, cex.lab=1.5 )

plot( delta.exo2, balance2, 
	xlab="delta estimate", ylab="Proportion of forward strand reads",
	cex.axis=1.5, cex.lab=1.5 )

plot( sigma.exo2, balance2, 
	xlab="sigma estimate", ylab="Proportion of forward strand reads",
	cex.axis=1.5, cex.lab=1.5 )

dev.off()

# plot

pdf( paste(dir.out,"delta_estimate_exo_seq_R1.pdf",sep="") )

### C 
par( mfrow=c(2,1) )
plot( density(delta.exo1), xlim=c(0,100),
	xlab="delta estimate", main="ChIP-exo",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
plot( density(delta.SET), xlim=c(0,100),
	xlab="delta estimate", main="ChIP-seq (SET)",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
dev.off()

pdf( paste(dir.out,"sigma_estimate_exo_seq_R1.pdf",sep="") )
### D
par( mfrow=c(2,1) )
plot( density(sigma.exo1), xlim=c(0,100),
	xlab="sigma estimate", main="ChIP-exo",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
plot( density(sigma.SET), xlim=c(0,100),
	xlab="sigma estimate", main="ChIP-seq (SET)",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
dev.off()

dr = "/p/keles/ChIPexo/volume3/ChIPexo/poster"

dt1 = data.table(delta = delta.exo1)
dt1[,Dataset :="ChIP-exo"]
dt2 = data.table(delta = delta.SET)
dt2[,Dataset:= "ChIP-seq (SET)"]
dt = rbind(dt1,dt2)
dt[ , Dataset := plyr::mapvalues( Dataset ,
        from = c("ChIP-exo","ChIP-seq (PET)","ChIP-seq (SET)"),
        to = c("exo","pet","set"))]

p1 = ggplot(dt,aes(delta,colour = Dataset))+geom_density(size=1)+facet_grid(Dataset~.)+
  theme_bw()+theme(legend.position = "none",plot.title = element_text(hjust = 0))+
  scale_color_manual( values = brewer.pal(3,"Set1")[c(1,3)])+
  ylab("")+xlab("delta estimate")
dt1 = data.table(sigma = sigma.exo1)
dt1[,Dataset :="ChIP-exo"]
dt2 = data.table(sigma = sigma.SET)
dt2[,Dataset:= "ChIP-seq (SET)"]
dt = rbind(dt1,dt2)

dt[ , Dataset := plyr::mapvalues( Dataset ,
        from = c("ChIP-exo","ChIP-seq (PET)","ChIP-seq (SET)"),
        to = c("exo","pet","set"))]

p2 = ggplot(dt,aes(sigma,colour = Dataset))+geom_density(size=1)+facet_grid(Dataset~.)+
  theme_bw()+theme(legend.position = "none",plot.title = element_text(hjust = 0))+
  scale_color_manual(values = brewer.pal(3,"Set1")[c(1,3)])+ylab("")+
  xlab("sigma estimate")

pdf(file = "figs/for_paper/sigma_delta_old_densities.pdf",width = 4 ,height  = 4)
print(p1 + ggtitle("C"))
print(p2 + ggtitle("D"))
dev.off()

pdf(file = file.path(dr,"sigma_estimate_exo_seq_R1.pdf"))
u = print(p1)
dev.off()


# strand balance

balance1 <- sapply( fit.exo1@fragSet, function(x) 
	if ( !is.na(x[1,1]) ) {
		return( mean( x[,3]=="F" ) )
	} else {
		return( NA )
	} )
balance1 <- balance1[ peak.R1[,5] > 3000 ]

balance2 <- sapply( fit.exo2@fragSet, function(x) 
	if ( !is.na(x[1,1]) ) {
		return( mean( x[,3]=="F" ) )
	} else {
		return( NA )
	} )
balance2 <- balance2[ peak.R2[,5] > 3000 ]

balanceSet <- sapply(fit.SET@fragSet ,function(x)
	if ( !is.na(x[1,1]) ) {
		return( mean( x[,3]=="F" ) )
	} else {
		return( NA )
	} )
nfrag = sapply(fit.SET@fragSet, nrow)
balanceSet = balanceSet[nfrag > 3e3]

dt1 = data.table(seq = "exo", fratio = balance1)
dt2 = data.table(seq = "set",fratio = balanceSet)
dt = rbind(dt1,dt2)

pdf(file = "figs/for_paper/forward_strand_ratio_comp_old.pdf",width = 4 ,height  = 4)
ggplot(dt , aes(fratio, colour = seq))+geom_density()+
  theme_bw()+scale_color_manual(values = brewer.pal(3,"Set1")[c(1,3)],name = "")+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Ratio of forward strand reads")+
  ylab("Density")+geom_vline(xintercept = .5,linetype = 2)+ggtitle("B")
dev.off()



pdf( paste(dir.out,"delta_estimate_exo_seq_R2.pdf",sep="") )
par( mfrow=c(2,1) )
plot( density(delta.exo2), xlim=c(0,100),
	xlab="delta estimate", main="ChIP-exo",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
plot( density(delta.SET), xlim=c(0,100),
	xlab="delta estimate", main="ChIP-seq (SET)",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
dev.off()

pdf( paste(dir.out,"sigma_estimate_exo_seq_R2.pdf",sep="") )
par( mfrow=c(2,1) )
plot( density(sigma.exo2), xlim=c(0,100),
	xlab="sigma estimate", main="ChIP-exo",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
plot( density(sigma.SET), xlim=c(0,100),
	xlab="sigma estimate", main="ChIP-seq (SET)",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
dev.off()


                     

# delta & sigma

delta.exo1 <- delta.exo1[ peak.R1[,5] > 3000 ]
sigma.exo1 <- sigma.exo1[ peak.R1[,5] > 3000 ]

delta.exo2 <- delta.exo2[ peak.R2[,5] > 3000 ]
sigma.exo2 <- sigma.exo2[ peak.R2[,5] > 3000 ]

delta.SET <- delta.SET[ !is.na(delta.SET) ]
sigma.SET <- sigma.SET[ !is.na(sigma.SET) ]

# plot (count >= 3k)

pdf( paste(dir.out,"count3k_issue_delta_estimate_exo_R1.pdf",sep="") )

a <- hexbin( delta.exo1, sigma.exo1, xbins = 30)
hbp <- plot( a, trans = log, inv = exp, 
	xlab = "delta estimate", ylab = "sigma estimate", 
	colramp = rainbow)

plot( delta.exo1, str1, log="y",
	xlab="delta estimate", ylab="Number of reads in the region",
	cex.axis=1.5, cex.lab=1.5 )

plot( sigma.exo1, str1, log="y",
	xlab="sigma estimate", ylab="Number of reads in the region",
	cex.axis=1.5, cex.lab=1.5 )

plot( delta.exo1, balance1, 
	xlab="delta estimate", ylab="Proportion of forward strand reads",
	cex.axis=1.5, cex.lab=1.5 )

plot( sigma.exo1, balance1, 
	xlab="sigma estimate", ylab="Proportion of forward strand reads",
	cex.axis=1.5, cex.lab=1.5 )

dev.off()

pdf( paste(dir.out,"count3k_issue_delta_estimate_exo_R2.pdf",sep="") )

a <- hexbin( delta.exo2, sigma.exo2, xbins = 30)
hbp <- plot( a, trans = log, inv = exp, 
	xlab = "delta estimate", ylab = "sigma estimate", 
	colramp = rainbow)

plot( delta.exo2, str2, log="y",
	xlab="delta estimate", ylab="Number of reads in the region",
	cex.axis=1.5, cex.lab=1.5 )

plot( sigma.exo2, str2, log="y",
	xlab="sigma estimate", ylab="Number of reads in the region",
	cex.axis=1.5, cex.lab=1.5 )

plot( delta.exo2, balance2, 
	xlab="delta estimate", ylab="Proportion of forward strand reads",
	cex.axis=1.5, cex.lab=1.5 )

plot( sigma.exo2, balance2, 
	xlab="sigma estimate", ylab="Proportion of forward strand reads",
	cex.axis=1.5, cex.lab=1.5 )

dev.off()

pdf( paste(dir.out,"count3k_delta_estimate_exo_seq_R1.pdf",sep="") )
par( mfrow=c(2,1) )
plot( density(delta.exo1), xlim=c(0,100),
	xlab="delta estimate", main="ChIP-exo",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
plot( density(delta.SET), xlim=c(0,100),
	xlab="delta estimate", main="ChIP-seq (SET)",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
dev.off()

pdf( paste(dir.out,"count3k_sigma_estimate_exo_seq_R1.pdf",sep="") )
par( mfrow=c(2,1) )
plot( density(sigma.exo1), xlim=c(0,100),
	xlab="sigma estimate", main="ChIP-exo",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
plot( density(sigma.SET), xlim=c(0,100),
	xlab="sigma estimate", main="ChIP-seq (SET)",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
dev.off()

pdf( paste(dir.out,"count3k_delta_estimate_exo_seq_R2.pdf",sep="") )
par( mfrow=c(2,1) )
plot( density(delta.exo2), xlim=c(0,100),
	xlab="delta estimate", main="ChIP-exo",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
plot( density(delta.SET), xlim=c(0,100),
	xlab="delta estimate", main="ChIP-seq (SET)",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
dev.off()

pdf( paste(dir.out,"count3k_sigma_estimate_exo_seq_R2.pdf",sep="") )
par( mfrow=c(2,1) )
plot( density(sigma.exo2), xlim=c(0,100),
	xlab="sigma estimate", main="ChIP-exo",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
plot( density(sigma.SET), xlim=c(0,100),
	xlab="sigma estimate", main="ChIP-seq (SET)",
	cex.axis=1.5, cex.lab=1.5, cex.main=1.5 )
dev.off()
