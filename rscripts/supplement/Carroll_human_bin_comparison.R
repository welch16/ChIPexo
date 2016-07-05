
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(dpeak)
library(viridis)
library(scales)
library(parallel)

devtools::load_all("~/Desktop/Docs/Code/ChIPUtils")


bin_size <- 150
frag_len <- 150

read_dr <- "/p/keles/ChIPexo/volume4/carroll_data/human/"

read_files <- list.files(read_dr,full.names = TRUE,recursive = TRUE,
                         pattern= "sort.bam")
read_files <- read_files[grep("bai",read_files,invert = TRUE)]

reads <- mcmapply(create_reads,read_files,FALSE, mc.cores = 6 )
names(reads) <- gsub(".sort.bam","",basename(read_files))

plots <- mapply(hexbin_plot,
  reads[1:3],reads[c(4,6,5)],
  MoreArgs = list(bin_size = bin_size, log = TRUE,frag_len = frag_len),
  SIMPLIFY = FALSE)#,mc.cores = 3)                  

plots <- lapply(plots,function(x){
  x + xlab("ChIP-Seq (SE) counts") + ylab("ChIP-exo counts")+
    theme_bw()+
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0))
})

plots <- mapply(function(x,y)x + ggtitle(y),plots,c("A","B","C"),SIMPLIFY = FALSE)
    
figs_dir <- "figs/supplement"
pdf(file = file.path(figs_dir,"carroll_ER_human_bin_comp.pdf"))
u <- lapply(plots,print)
dev.off()
