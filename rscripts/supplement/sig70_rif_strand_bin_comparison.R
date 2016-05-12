
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(dpeak)
library(ChIPUtils)
library(viridis)
library(scales)
library(parallel)

bin_size <- 150
frag_len <- 150

read_dr <- "/p/keles/ChIPexo/volume7/Landick/K12"

read_files <- list.files(read_dr,full.names = TRUE,recursive = TRUE,
                         pattern= "sort.bam")
read_files <- read_files[grep("bai",read_files,invert = TRUE)]
read_files <- read_files[grep("seed",read_files,invert = TRUE)]
read_files <- read_files[grep("rif",read_files)]

read_files <- read_files[grep("SET",read_files,invert = TRUE)]
read_files <- read_files[grep("1369",read_files,invert = TRUE)]

## reads <- mclapply(read_files,readGAlignments,param = NULL,mc.cores = 8)
## reads <- mclapply(reads,as,"GRanges",mc.cores = 8)
## names(reads) <- gsub(".sort.bam","",basename(read_files))

library(devtools)
load_all("~/Desktop/Docs/Code/ChIPUtils")

reads <- mcmapply(create_reads,read_files,rep(c(FALSE,TRUE),each = 4),
                  mc.cores = 8 )

plots <- mcmapply(hexbin_plot,
  reads[5:8],reads[1:4],
  MoreArgs = list(bin_size = bin_size, log = TRUE,frag_len = frag_len),
  SIMPLIFY = FALSE,mc.cores = 4)                  

plots <- lapply(plots,function(x){
  x + xlab("ChIP-Seq (PE) counts") + ylab("ChIP-exo counts")+
    theme_bw()+
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0))
})

plots <- mapply(function(x,y)x + ggtitle(y),plots,c("A","B","C","D"),SIMPLIFY = FALSE)
    
figs_dir <- "figs/supplement"
pdf(file = file.path(figs_dir,"sig70_rif_count_comparison.pdf"))
u <- lapply(plots,print)
dev.off()
