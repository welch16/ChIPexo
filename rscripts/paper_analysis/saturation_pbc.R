
rm(list = ls())

library(ChIPUtils)
library(GenomicAlignments)
library(data.table)
library(parallel)

## pet pbc is really high, even more given that in E.Coli, the pbc is usually smaller than eukaryote genomes

indir <- "/p/keles/ChIPexo/volume6/K12/saturation"
files <- list.files(indir,recursive = TRUE)
files <- files[grep("samples",files)]
files <- files[grep("bai",files,invert = TRUE)]

exo <- files[grep("exo",files)]
pet <- files[grep("PET",files)]
set <- files[grep("SET",files)]

get_pbc <- function(file,is_PET = FALSE){

  reads <- create_reads(file,is_PET)
  PBC(reads)
}

exo_pbc <- mclapply(file.path(indir,exo),get_pbc,is_PET = FALSE,mc.cores = detectCores())
names(exo_pbc) <- basename(exo)

pet_pbc <- mclapply(file.path(indir,pet),get_pbc,is_PET = TRUE,mc.cores = detectCores())
names(pet_pbc) <- basename(pet)

set_pbc <- mclapply(file.path(indir,set),get_pbc,is_PET = FALSE,mc.cores = detectCores())
names(set_pbc) <- basename(set)

convert_to_table <- function(pbc)
{
  info <- names(pbc)
  info <- strsplit(info,"_",fixed = TRUE)
  out <- data.table(edsn = sapply(info,function(x)x[1]),
                    samp = sapply(info,function(x)x[3]),
                    seed = sapply(info,function(x)x[4]))
  out[,edsn := gsub("edsn","",edsn)]
  out[,samp := gsub("samp","",samp)]
  out[,samp := as.numeric(gsub("K","",samp))*1e3]
  out[,seed := gsub(".bam","",seed)]
  out[,pbc := do.call(c,pbc)]
  return(out)                
}

exo_pbc <- convert_to_table(exo_pbc)
pet_pbc <- convert_to_table(pet_pbc)
set_pbc <- convert_to_table(set_pbc)

pbc <- list(exo = exo_pbc,pet = pet_pbc,set = set_pbc)

save(pbc, file = "data/saturation_pbc.RData")

load("data/saturation_pbc.RData")
pbc <- mapply(function(x,y)x[,protocol := y],pbc,names(pbc),SIMPLIFY = FALSE)
pbc <- do.call(rbind,pbc)
pbc[,rif := ifelse(edsn %in% c(1311,1317,1396,1400),"0min","20min")]
pbc[,repl := ifelse(edsn %in% c(1311,1314,1396,1398),"Rep-1","Rep-2")]


pdf(file = "figs/saturation/K12_alignment/PBC_saturation.pdf")
ggplot(pbc[,mean(pbc),by = .(samp,protocol,rif,repl)],aes(samp,V1,colour = protocol))+geom_line(size = 1.5)+
  facet_grid(repl ~ rif)+theme_bw()+theme(legend.position = "top")+
  scale_color_brewer(palette = "Set1")+xlab("Nr. of reads")+ylab("Average PBC")
dev.off()
  
