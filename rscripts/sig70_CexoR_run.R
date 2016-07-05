
rm(list = ls())

library(CexoR)
library(data.table)
library(GenomicAlignments)

## file <- "/p/keles/ChIPexo/volume6/sigma70_old_alignment/ChIP-exo_sigma70_exp_phase_R1.sort.bam"

file <- "/p/keles/ChIPexo/volume7/Landick/K12/ChIPexo/aerobic_vs_anaerobic/edsn931_Sig70.sort.bam"

gr <- readGAlignments(file,param = NULL)
nms <- unique(as.character(seqnames(gr)))

interval <- GRanges(seqnames = nms , range = IRanges(1,1e6))
gr2 <- readGAlignments(file , param = ScanBamParam(which  = interval))

size <- read.table("/p/keles/ChIPexo/volume7/Landick/K12/K12_size")
peak_pair <- cexor(bam = file, chrN = nms,
    chrL = size$V2,bedfile = FALSE,
    N = 4e6,p = 1e-2)
