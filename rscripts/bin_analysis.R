


rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(dpeak)
library(ChIPUtils)
library(viridis)
library(scales)


## devtools::load_all("~/Desktop/Docs/Code/ChIPUtils")


######################################################################################

## Condition table

source("R/base_edsn.R")

what <- c("exo","pet","set")
char <- lapply(what,edsn_tab_old)
names(char) <- what


#######################################################################################

### reads directories

chip_dirs <- list()
chip_dirs[["exo"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"
chip_dirs[["pet"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET"
chip_dirs[["set"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET"

######################################################################################

## get files

read_files <- mapply(filter_reads,
  chip_dirs,char,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)

reads <- mapply(create_reads,read_files,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)


pdf(file = file.path("figs/for_paper","ChIPseqPET_ChIPexo_tagCount_comparison.pdf"),width = 4,height = 4)
p <- hexbin_plot(reads[[2]],reads[[1]],150,frag_len = 150)+xlab("ChIP-seq (PE) counts")+ylab("ChIP-exo counts")+theme_bw()
dt <- copy(p$data)
dt[,x := 1 + x]
dt[,y := 1 + y]
print(p)    
print(p + xlim(0,500)+ylim(0,500)+coord_fixed())
z <- p %+% dt +coord_fixed() + scale_x_log10() + scale_y_log10()+
  theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  ggtitle("A")
print(z)
dev.off()

reg <- GRanges(seqnames = "U00096",ranges = IRanges(start = 1, end = 4641652))
bins <- create_bins(150,chrom = reg)

count_strand <- function(reads,bins)
{
  fwd <- unlist(lapply(readsF(reads),dt2gr))[[1]]
  bwd <- unlist(lapply(readsR(reads),dt2gr))[[1]]
  f <- countOverlaps(bins,fwd)
  r <- countOverlaps(bins,bwd)
  out <- data.table(f,r)
  out <- cbind(gr2dt(bins),out)
  out <- out[ f > 0 & r > 0]
  out[, M := log2(f * r)]
  out[, A := log2(f /r)]
  out[,ratio:= f / (f + r)]
  return(out)
}

strand_count <- lapply(reads,count_strand,bins)
strand_count <- mapply(function(x,u)x[,what := u],strand_count,c("exo","pet","set"),SIMPLIFY = FALSE)
strand_count <- do.call(rbind,strand_count)

filter <- 100

strand_count_filter <- copy(strand_count)
strand_count_filter <- strand_count_filter[f + r > filter]

strand_count[ , which := "all bins"]

strand_count_filter[ , which := paste0("nr. reads > ",filter)]
out <- rbind(strand_count,strand_count_filter)


r <- viridis(100,option = "D")

pdf(file = file.path("figs/for_paper","MA_plot_fwd_bwd_by_seq.pdf"),width = 8,height = 5)
ma <- ggplot(strand_count,aes(M , A))+stat_binhex(bins = 100)+
  scale_fill_gradientn(colours = r,trans = "log10",labels = trans_format("log10",math_format(10^.x)))+
  facet_grid(which ~ what)+theme_bw()+theme(legend.position = "top")+coord_fixed(xlim = c(0,30), ylim = c(-12,12))
print(ma)
print(ma %+% strand_count_filter)
print(ma %+% out)     
dev.off()      


plots <- list(p,z,ma)

save(plots,file = "data/for_paper/ChIPseqPET_ChIPexo_tagCount_comparison.RData")




