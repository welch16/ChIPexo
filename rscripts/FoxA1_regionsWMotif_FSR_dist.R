
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)

library(ggplot2)
library(RColorBrewer)

library(magrittr)
library(dplyr)

## load_all("~/Desktop/Docs/Code/Segvis")
## load_all("~/Desktop/Docs/Code/ChIPexoQual")

data_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"
data_dir <- "data/ChIPexo_QC_runs"

files <- list.files(data_dir,full.names = TRUE,include.dirs = TRUE)
files <- files[grep("mouse",files)]

load_file <- function(x){
  load(x)
  return(ext_stats[["stats"]])
}

stats <- mclapply(files,load_file,mc.cores = 3)
names(stats) <- gsub(".RData","",basename(files))

## exo <- lapply(files,create_exo_experiment,parallel = TRUE)
## stats <- lapply(exo,summary_stats)

## exo_reads <- lapply(exo,reads)
## readlength <- mclapply(exo_reads,function(x)table(width(x)),
##                 mc.cores =3)

rl <- 36 ## most repeated read length

### param to get a high quality set of regions

## 1 - only regions with reads in both strands
stats <- lapply(stats,function(x)x[f > 0 & r > 0])

## 2 - wider regions
stats <- lapply(stats,function(x)x[between(width,3 * rl,1e3)])
                                   
## 3 - deeper regions
minNpos <- 15
minDepth <- 100

stats <- lapply(stats,function(x)x[npos > minNpos])
stats <- lapply(stats,function(x)x[depth > minDepth])

## remove chrM regions
stats <- lapply(stats,function(x)x[seqnames != "chrM"])




exo_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo"

fimo_files <- list.files(exo_dir,pattern = "txt",recursive = TRUE,
                         include.dirs = TRUE,full.names = TRUE)

fimo <- mclapply(fimo_files,read.table,mc.cores = 12)
fimo <- lapply(fimo,data.table)
names(fimo) <- basename(dirname(fimo_files))

##     both only_reg only_peak  rep
## 1: 12726      991       502 rep1
## 2:  3271      332       442 rep2
## 3:  4074       32      3049 rep3

# V1 motif id
# V2 sequence id
# V3 start
# V4 end
# V5 strand
# V6 score
# V7 p.value
# V8 q.value
# V9 sequence

fimo <- lapply(fimo,function(x){
  setnames(x,names(x),c("motifID","sequenceID",
                        "motifStart","motifEnd","strand",
                        "score","pval","qval","sequence"))
  return(x)})

fimo <- fimo[grep("FOXA1",names(fimo))]

chr_from_fimo <- function(x)sapply(strsplit(as.character(x),":"),
                function(y)y[1])

start_from_fimo <- function(x){
  out <- sapply(strsplit(as.character(x),":"),function(y)y[2])
  out <- sapply(strsplit(out,"-",fixed = TRUE),function(z)z[1])
  return(as.numeric(out))
}

end_from_fimo <- function(x){
  out <- sapply(strsplit(as.character(x),":"),function(y)y[2])
  out <- sapply(strsplit(out,"-",fixed = TRUE),function(z)z[2])
  return(as.numeric(out))
}

fimo <- lapply(fimo,function(x){
  x[,seqnames := chr_from_fimo(sequenceID)]
  x[,start := start_from_fimo(sequenceID)]
  x[,end := end_from_fimo(sequenceID)]
  return(x)})

chr2GR <- function(cc)
{
    cc1 = cc %>% strsplit("\\:")
    cc2 = cc1 %>% sapply(function(x)x[2])
    
    GRanges(
        seqnames = cc1  %>% sapply(function(x)x[1]),
        ranges =
            IRanges(
                start = cc2 %>% strsplit("\\-") %>% sapply(function(x)x[1]) %>% as.numeric,
                end = cc2 %>% strsplit("\\-") %>% sapply(function(x)x[2]) %>% as.numeric))
    

}

match_by_name <- function(fimo,stats)
{
    c1 = fimo$sequenceID %>% as.character 
    c2 = stats$match %>% as.character
    c1 = gsub("-","_",c1)
    setkey(stats,match)

    stats[c1] %>% as.tbl
        
}

stats_fimo = list()
stats_fimo[["rep1"]] = match_by_name(fimo[[2]],stats[[1]])
stats_fimo[["rep2"]] = match_by_name(fimo[[3]],stats[[2]])
stats_fimo[["rep3"]] = match_by_name(fimo[[1]],stats[[3]])


stats_fimo_summary <- function(cc)
{
    cc %>% group_by(match) %>%
        summarize(
            nmotif = n(),
            fsr = mean(fsr),
            depth = mean(depth),
            npos = mean(npos),
            width = mean(width)) %>% ungroup
}

summaries = lapply(stats_fimo,stats_fimo_summary)

theme_set(theme_bw())


figs = "figs/NAR_review/FSR_peaks"

pdf(file.path(figs,"FSR_proportion_plot_regionsWMotif.pdf"))
plots = lapply(summaries,
       function(x){
           x %>%
               ggplot(aes(fsr,fill = as.factor(nmotif)))+
                 geom_histogram(aes(y = ..count.. / sum(..count..)),bins = 35,colour = "black" )+
                 scale_fill_brewer(palette = "Pastel2",name = "# motifs in region")+
                 theme(axis.title.x = element_blank(),legend.position = "top")+
                 ylab("Proportion of regions")+
                     geom_vline(xintercept = .5,linetype = 2)})
u = lapply(plots,print)
dev.off()

table = lapply(summaries,
               function(x){
                   c("nmotif" = x %>% select(nmotif) %>% sum,
                     "fsr_le_0.1" = x %>% filter(fsr <= 0.1) %>% select(nmotif) %>% sum,                     
                     "fsr_le_0.25" = x %>% filter(fsr <= 0.25) %>% select(nmotif) %>% sum,
                     "fsr_ge_0.75" = x %>% filter(fsr >= 0.75) %>% select(nmotif) %>% sum,
                     "fsr_ge_0.9" = x %>% filter(fsr >= 0.9) %>% select(nmotif) %>% sum)
                   })
table = do.call(rbind,table) %>%
    as.data.frame %>% as.tbl %>% mutate(Rep = paste0("Rep-",seq_len(3))) %>%
    select(Rep,everything())

table2 = table %>%
    mutate(fsr_le_0.1 = fsr_le_0.1 / nmotif * 100,
           fsr_le_0.25 = fsr_le_0.25 / nmotif * 100,
           fsr_ge_0.75 = fsr_ge_0.75 / nmotif * 100,
           fsr_ge_0.9 = fsr_ge_0.9 / nmotif * 100,
           nmotif = NULL)
