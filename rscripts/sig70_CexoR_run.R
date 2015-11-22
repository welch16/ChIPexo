
rm(list = ls())

library(CexoR)
library(data.table)
library(GenomicAlignments)

file <- "/p/keles/ChIPexo/volume6/sigma70_old_alignment/ChIP-exo_sigma70_exp_phase_R1.sort.bam"

gr <- readGAlignments(file,param = NULL)
nms <- unique(as.character(seqnames(gr)))

interval <- GRanges(seqnames = nms , range = IRanges(1,1e6))
gr2 <- readGAlignments(file , param = ScanBamParam(which  = interval))


peak_pair <- cexor(bam = file, chrN = nms,chrL = 1e7,bedfile = FALSE,dpeaks = c(0,200),dpairs = 200,
   N = 1e3,p = 1)
