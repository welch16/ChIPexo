


rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(dpeak)
library(ChIPUtils)

## devtools::load_all("~/Desktop/Docs/Code/ChIPUtils")


######################################################################################

## Condition table

source("R/base_edsn.R")

what <- c("exo","pet","set")
char <- lapply(what,edsn_tab_old)
names(char) <- what

bin_size

#############1#########################################################################

### reads directories

chip_dirs <- list()
chip_dirs[["exo"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"
chip_dirs[["pet"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_PET"
chip_dirs[["set"]] <- "/p/keles/ChIPexo/volume3/LandickData/ChIPseq_SET"

######################################################################################

## get files

filter_reads <- function(dir,char,isPET)
{
  files <- list.files(dir)
  files <- files[grep(char[,(edsn)],files)]
  files <- files[grep("bai",files,invert = TRUE)]  
  if(isPET){
    files <- files[grep("filter",files)]
  }
  files <- files[grep("run",files,invert = TRUE)]
  return(file.path(dir,files))
}

read_files <- mapply(filter_reads,
  chip_dirs,char,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)

reads <- mapply(create_reads,read_files,c(FALSE,TRUE,FALSE),SIMPLIFY = FALSE)

pdf(file = file.path("figs/for_paper","ChIPseqPET_ChIPexo_tagCount_comparison.pdf"))
p <- hexbin_plot(reads[[2]],reads[[1]],150,frag_len = 150)+xlab("ChIP-seq (PET) tag counts")+ylab("ChIP-exo tag counts")+theme_bw()
dt <- copy(p$data)
dt[,x := 1 + x]
dt[,y := 1 + y]
p    
p + xlim(0,500)+ylim(0,500)+coord_fixed()
p %+% dt +coord_fixed() + scale_x_log10() + scale_y_log10()
dev.off()

