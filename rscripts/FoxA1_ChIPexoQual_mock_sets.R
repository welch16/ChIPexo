
rm(list = ls())

library(GenomicAlignments)
library(GenomicRanges)
library(data.table)

indir = "/p/keles/ChIPexo/volume4/carroll_data/mouse"
files = list.files(indir,full.names = TRUE,pattern = "sort")
files = files[grep("bai",files,invert = TRUE)]

outdr = indir
outfiles = file.path(outdr,
  paste0("ChIPexo_carroll_FoxA1_mouse_rep", c(3,1,2),"_first3chr.bam"))
outfiles1 = file.path(outdr,
  paste0("ChIPexo_carroll_FoxA1_mouse_rep",c(3,1,2),"_chr1.bam"))



filter_factory <- function(){
  list()
}

sizes = fread("/p/keles/SOFTWARE/hg19.chrom.sizes")

gr = sizes[1:3,GRanges(seqnames = V1,
  range = IRanges(start = 1,end = V2))]

gr1 = sizes[1,GRanges(seqnames = V1,
  range = IRanges(start = 1,end = V2))]


a <- mapply(function(infi,outfi){
  filter = FilterRules(filter_factory())
  filterBam(infi,outfi,filter = filter,
            param = ScanBamParam(which = gr))
},files,outfiles,SIMPLIFY =  FALSE)

a1 <- mapply(function(infi,outfi){
  filter = FilterRules(filter_factory())
  filterBam(infi,outfi,filter = filter,
            param = ScanBamParam(which = gr1))
},files,outfiles1,SIMPLIFY =  FALSE)
