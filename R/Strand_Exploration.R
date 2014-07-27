
rm(list = ls())
library(devtools)
library(GenomicAlignments)
library(parallel)



dr ="/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo3/rawdata"
folder = c("ChIPexo","ChIPseq_PET","ChIPseq_SET")

## Set bam files on list
files = lapply(folder,function(x,dr){
  ff = list.files(file.path(dr,x))
  return(ff[!grepl("bai",ff) & !grepl("sam",ff)])},dr)
names(files) = folder

ff = lapply(folder,function(x,dr,files)file.path(dr,x,files[[x]]),dr,files)
names(ff) = folder

fragLen = 0

strand_coverage <- function(reads,peak,st, fragLen = 0)
{
  reads = as(reads,"GRanges")
  reads = subset(reads, subset = strand(reads) == st)
  if(fragLen >0){
    reads = resize(reads,fragLen)
  }  
  match = subjectHits(findOverlaps(peak,reads))
  cover = coverage(reads[match])[[1]]
  x = seq(start(peak),end(peak))
  if(nrun(cover) == 1){
    y = rep(runValue(cover),length(x))
  }else{
    xp = cumsum(runLength(cover)[1:(nrun(cover)-1)])
    yp = runValue(cover)
    y = stepfun(xp,yp)(x)
  }
  return(y)  
}

mainFromPeak <- function(peak)paste0(seqnames(peak),":",start(peak),"-",end(peak))

plot_strand <- function(reads,peak,main = "",fragLen=0)
{
  x = seq(start(peak),end(peak))
  yF = suppressWarnings(strand_coverage(reads,peak,"+"))
  yR = suppressWarnings(strand_coverage(reads,peak,"-"))
  yl = c(0, max(max(yF),max(yR))*1.2)

  par( mar=c(4,4,0.5,0.5), oma=c(0,0,2,0),mfrow = c(3,1) )   
  plot(x,yF,type = "l",col = "red",,ylim = c(0,max(yF)*1.2) , xlab = "Genomic coordinates F-strand",ylab = "Counts")  
  plot(x,yR,type = "l",col = "blue",ylim = c(0,max(yR)*1.2),xlab = "Genomic coordinates R-strand",ylab = "Counts" )
  plot(x,yF,type = "l",col = "red",,ylim = yl , xlab = "Genomic coordinates",ylab = "Counts")
  lines(x,yR,col = "blue")
  mtext( main,outer = TRUE)
  
}

#pdf("test.pdf")
#lapply(set,FUN = plot_strand,peak,main = "ADGD")
#dev.off()

## start = c(576875,798921,1152674,1192831,1820005,2310506,2996662,3351853,3973436,4262028)
## end = c(577492,799434,1153291,1193448,1820696,2311153,2997175,3352366,3973949,4262541)
## peaks = GRanges(seqnames = "U00096",ranges = IRanges(start = start,end = end),strand = "*")
## width(peaks) = 1000

exo = mclapply(ff[[1]],FUN = readGAlignmentsFromBam,param = NULL,mc.cores = 8)  

set = mclapply(ff[[3]],FUN = readGAlignmentsFromBam,param = NULL,mc.cores = 8)

pet_flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isUnmappedQuery = NA,
                hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                isFirstMateRead = NA, isSecondMateRead = NA, isNotPrimaryRead = NA,
                isNotPassingQualityControls = NA, isDuplicate = NA)
pet_param = ScanBamParam(flag = pet_flag, simpleCigar = FALSE,
                 reverseComplement = FALSE, tag = character(0),
                 what = character(0))
pet = mclapply(ff[[2]][-c(11,12)],FUN = readGAlignmentsFromBam,param = pet_param,mc.cores = 8)

save(list = c("exo","pet","set"),file = "../RData/Strand.RData")

fileFromPeak <- function(seq,peak) paste(seq,as.character(seqnames(peak)),start(peak),end(peak),sep = "_") 

names(exo) = sub("_qc.sorted.bam","",files[[1]])
names(pet) = sub("_qc.sorted.bam","",files[[2]][-(11:12)])
names(set) = sub("_qc.sorted.bam","",files[[3]])

## N = length(peaks)
## figsDr = "../Figs/Strand"
## mclapply(1:N,function(i,peaks,exo,pet,set){
##   pdf(file = file.path(figsDr,paste0(fileFromPeak("ChIP-Exo",peaks[i]),".pdf")))
##      lapply(1:length(exo),function(j,exo,peak)plot_strand(exo[[j]],peak,main = names(exo)[j]),exo,peaks[i])
##   dev.off()
##   pdf(file = file.path(figsDr,paste0(fileFromPeak("ChIP-Seq-PET",peaks[i]),".pdf")))
##     lapply(1:length(pet),function(j,pet,peak)plot_strand(pet[[j]],peak,main = names(pet)[j]),pet,peaks[i])
##   dev.off()
##   pdf(file = file.path(figsDr,paste0(fileFromPeak("ChIP-Seq-SET",peaks[i]),".pdf")))
##     lapply(1:length(set),function(j,set,peak)plot_strand(set[[j]],peak,main = names(set)[j]),set,peaks[i])
##   dev.off()
## },peaks,exo,pet,set,mc.cores = 8)

library(gdata)
tab = list()
tab[[1]] = read.xls("../Alignment/ChIP-Exo-summary.xls",sheet = 1)[,1:7]
tab[[2]] = read.xls("../Alignment/ChIP-Exo-summary.xls",sheet = 2)[,1:7]
tab[[3]] = read.xls("../Alignment/ChIP-Exo-summary.xls",sheet = 3)[,1:7]
tab[[1]]$seq = "Exo"
tab[[2]]$seq = "PET"
tab[[3]]$seq = "SET"
nn = names(tab[[1]])
tab = do.call(rbind,tab)
tab = do.call(cbind,lapply(1:ncol(tab),function(i,tab){
  if(class(tab[,i]) == "factor"){
    x = as.character(tab[,i])
    x[x == ""] = NA
  }else{
    x = tab[,i]
  }
return(x)},tab))
tab = as.data.frame(tab)
names(tab) = nn
rownames(tab)=NULL

resume.samples <- function(edsn = NULL,cult = NULL,ip = NULL,phase = NULL,growth = NULL,
  rif = NULL,rep = NULL,seq = NULL)
{
  empty = TRUE
  st = ""
  if(!is.null(edsn)){st = paste0(st,ifelse(!empty,"&",""),"edsn=='",edsn,"'");empty = FALSE}
  if(!is.null(cult)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(cult),"is.na(cult)",paste0("cult==",cult)));empty = FALSE}
  if(!is.null(ip)){st = paste0(st,ifelse(!empty,"&",""),"ip=='",ip,"'");empty = FALSE}
  if(!is.null(phase)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(phase),"is.na(phase)",paste0("phase=='",phase,"'")))
    empty = FALSE}
  if(!is.null(growth)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(growth),"is.na(growth)",paste0("growth=='",growth,"'")))
    empty = FALSE}  
  if(!is.null(rif)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(rif),"is.na(rif)",paste0("rif=='",rif,"'")));empty = FALSE}
  if(!is.null(rep)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(rep),"is.na(rep)",paste0("rep==",rep)));empty = FALSE}
  if(!is.null(seq))st = paste0(st,ifelse(!empty,"&",""),"seq=='",seq,"'") 
  return(st)
}

bin.density <- function(bins,reads)
{  
  counts_F = countOverlaps(bins,subset(reads,subset = strand(reads) == "+"))
  counts_R = countOverlaps(bins,subset(reads,subset = strand(reads) == "-"))
  ratio = (counts_F + 1)/(counts_F + counts_R + 2) # used pseudo-counts to avoid zero denominator
  return(density(ratio))
}


plot.density <- function(binSize,exo.sets,pet.sets,genomeLength = seqlengths(exo.sets[[1]]),main = "")
{
  bins = GRanges(seqnames = names(genomeLength),ranges = IRanges(start = seq(1,genomeLength,by=binSize),width = binSize),strand = "*")
  exo.densities = suppressWarnings(mclapply(exo.sets,function(x,bins)bin.density(bins,x),bins,mc.cores = 2))
  pet.densities = suppressWarnings(mclapply(pet.sets,function(x,bins)bin.density(bins,x),bins,mc.cores =2))
  par( mar=c(4,4,0.5,0.5), oma=c(0,0,2,0),mfrow = c(2,1) )   
  yl = c(0,max(exo.densities[[1]]$y,exo.densities[[2]]$y,pet.densities[[1]]$y,pet.densities[[2]]$y)*1.2)
  xl = c(0,1)
  plot(exo.densities[[1]]$x ,exo.densities[[1]]$y,type = "l",ylim = yl,xlim = xl,xlab = "Fwd. strand ratio",ylab = "Density")
  lines(pet.densities[[1]]$x , pet.densities[[1]]$y,col = "red")
  legend("topright",lty=c(1,1),c("ChIP-Exo","ChIP-Seq-PET"),col=c("black","red")) 
  plot(exo.densities[[2]]$x,exo.densities[[2]]$y,type = "l",ylim = yl,xlim = xl,xlab = "Fwd. strand ratio",ylab = "Density")
  lines(pet.densities[[2]]$x,pet.densities[[2]]$y,col = "red")
  legend("topright",lty=c(1,1),c("ChIP-Exo","ChIP-Seq-PET"),col=c("black","red"))   
  mtext( main,outer = TRUE)

}


ip = c("Sig70","BetaPrimeFlag")
rif = c("0 min","20 min")
growth = "Aerobic"
phase = "Exponential"
j = 1
st = list()
for(i in ip){
  for(r in rif){
    st[[j]] = resume.samples(ip = i,rif = r,growth = growth,phase = phase)
    j=j+1
  }
}
    
binSize = c(25,50,100,200,500,100)

for(bin in binSize){
  pdf(file = file.path("../Figs/Densities/",paste0("Densities_RifExperiment_",bin,".pdf")))
  lapply(st1,function(x,bin,tab,exo,pet){
    tt = subset(tab,subset = eval(parse(text =x)))
    edsn = as.character(tt$edsn)
    exo.sets = names(exo)[do.call(c,lapply(edsn,FUN = grep,names(exo)))]
    exo.sets = lapply(exo.sets,function(y,exo)exo[[y]],exo)
    pet.sets = names(pet)[do.call(c,lapply(edsn,FUN = grep,names(pet)))]
    pet.sets = lapply(pet.sets,function(y,pet)pet[[y]],pet)
    plot.density(bin,exo.sets,pet.sets,main = x)
  },bin,tab,exo,pet)
  dev.off()    
}

                         




