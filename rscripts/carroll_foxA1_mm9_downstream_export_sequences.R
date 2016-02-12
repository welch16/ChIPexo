
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(BSgenome.Mmusculus.UCSC.mm9)

sitedir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/sites"
files <- list.files(sitedir)

sites <- lapply(file.path(sitedir,files),read.table,header = TRUE)
sites <- lapply(sites,data.table)

bandwidth <- 50

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
