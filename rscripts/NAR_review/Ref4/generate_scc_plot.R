#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--files",action = "store_true",type = "character",
                help = "Aligned read files to generate SCC curves"),
    make_option("--names",action = "store_true",type = "character",
                help = "Names of the aligned reads samples"),
    make_option("--figs",action = "store_true",type = "character",default = "./Rplots",
                help = "Figures directory + prefix")
    )

opt = parse_args(OptionParser(option_list = optList))

library(ChIPQC)
library(parallel)

separateFiles <- function(ff)
{
  if(grepl(",",ff)){
      ff = ff %>% strsplit(',') %>% unlist     
  }else{
      ff = Sys.glob(ff)      
  }
  ff
}

opt$files = separateFiles(opt$files)

chipqc = mclapply(opt$files,ChIPQCsample,annotation = NULL,runCrossCor = TRUE,mc.cores = 10)
opt$names = opt$names %>% strsplit(",") %>% unlist

merge_scc <- function(chipqc,names)
{
    scc = lapply(chipqc , function(x)x@CrossCorrelation)
    scc = mapply(function(x,y){
        x = tibble(cross.cor = x)
        x %>% mutate(
                  shift = seq_len(nrow(x)),
                  sample = y) %>% select(sample,shift,cross.cor)},
        scc,names,SIMPLIFY = FALSE)
    bind_rows(scc)                                                          
}


scc = merge_scc(chipqc,opt$names)

theme_set(theme_bw())

pdf(paste0(opt$figs,"_scc.pdf"),height = 7,width = 9)
p = ggplot(scc ,aes(shift,cross.cor,colour = sample)) + geom_line() +
    scale_color_brewer(palette = "Set1") +
    theme(legend.position = "top")
print(p)
dev.off()
