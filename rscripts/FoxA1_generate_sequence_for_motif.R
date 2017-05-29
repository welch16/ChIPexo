
rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(devtools)
library(parallel)

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

## convert to regions

regions <- lapply(stats,function(x)
                  ChIPUtils::dt2gr(x[,2:4,with = FALSE]))
regions <- lapply(regions,sort)

outdir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse"

## check that the "high quality" regions overlap with peaks
peaks <- lapply(list.files(outdir,
   pattern = "peaks",recursive = TRUE,full.names= TRUE),
               read.table)
peaks <- lapply(peaks,data.table)
peaks <- lapply(peaks,function(x)x[,1:3,with = FALSE])
peaks <- lapply(peaks,function(x){
  setnames(x,names(x),c("seqnames","start","end"))
  return(ChIPUtils::dt2gr(x))})

venn <- mapply(function(region,peak){
  ov <- findOverlaps(region,peak)
  inter_peak <- length(unique(subjectHits(ov)))
  inter_reg <- length(unique(queryHits(ov)))
  out <- data.table(both = inter_reg,
                    only_reg = length(region) - inter_reg,
                    only_peak = length(peak) - inter_peak)
  return(out)                   
},regions,peaks,SIMPLIFY = FALSE)
venn <- do.call(rbind,venn)

venn[, rep := c("rep3","rep1","rep2")]

venn <- venn[order(rep)]

##     both only_reg only_peak  rep
## 1: 12726      991       502 rep1
## 2:  3271      332       442 rep2
## 3:  4074       32      3049 rep3


## must of the high quality regions are peaks, so we are getting
## high quality regions (i.e. good sanity check)

library(BSgenome.Mmusculus.UCSC.mm9)

sequences <- mclapply(regions,function(x)
   getSeq(Mmusculus,x),mc.cores = 3)

sequences <- mclapply(sequences,function(x)as.character(x),mc.cores =3)

nms <- lapply(regions,function(x){
  paste0(as.character(seqnames(x)),":",start(x),"-",end(x))})

fasta_formats <- mapply(function(nms,seqs)paste0(">",nms,"\n",seqs),
                        nms,sequences,SIMPLIFY = FALSE)

exo_dir <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/fasta_format_exo"

mapply(write.table,fasta_formats,
       file.path(exo_dir,"sequences",
                 gsub("sort.bam","sequences.fna",basename(files))),
  MoreArgs = list(quote = FALSE,row.names = FALSE,col.names = FALSE))




## library(rGADEM)
## library(MotIV)

## data(FOXA1_rGADEM)
## motifs <- getPWM(gadem)

## runs <- data.table(expand.grid(1:length(regions),1:length(motifs)))

## gademList <- mcmapply(function(i,j,regions,motifs){
##   region <- regions[[i]]
##   motif <- motifs[[j]]
##   gad <- GADEM(as(region,"RangedData"),
##                genome = Mmusculus, seed = motif)
##   return(gad)},runs[,(Var1)],runs[,(Var2)],
##                     MoreArgs = list(regions,motifs),SIMPLIFY = FALSE,
##                     mc.cores = 21)
                   
## save(gademList , file = "data/FoxA1_ChIPexo_rGADEM_analysis.RData")                    

