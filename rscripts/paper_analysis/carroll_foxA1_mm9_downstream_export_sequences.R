
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(BSgenome.Mmusculus.UCSC.mm9)

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

bandwidth <- 10

anchors <- lapply(sites,function(x,bw){
  summits <- x[,mid(IRanges(start = start,end = end))]
  out <- GRanges(seqnames = x[,(chrID)],
                 ranges = IRanges(
                   start = summits - bw ,
                   end = summits + bw))
  return(out)
},bandwidth)

sequences <- lapply(anchors,function(x)
   getSeq(Mmusculus,x,as.character = TRUE))

nms <- lapply(sites,function(x)x[,paste0(chrID,":",start,"-",end)])


fasta_formats <- mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),nms,sequences,SIMPLIFY = FALSE)

mapply(write.table,fasta_formats,file.path(dirname(sitedir),"sequences",gsub("sites.txt","sequences.fna",files)),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))
