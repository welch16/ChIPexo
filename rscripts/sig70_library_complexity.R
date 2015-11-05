
rm(list = ls())

library(GenomicAlignments)
library(data.table)
library(ChIPUtils)
library(ggplot2)

exo_dir <- "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"

files <- list.files(exo_dir)

files <- files[grep("bam",files)]
files <- files[grep("edsn",files)]
files <- files[grep("bai",files,invert = TRUE)]

source("R/base_edsn.R")

exo <- list(edsn_tab("exo"),edsn_tab_old("exo"))
exo[[2]][,growth := NULL]
exo[[3]] <- copy(exo[[2]])
exo[[3]][,edsn := 933]
exo[[3]][,repl := 2]
exo[[1]] <- exo[[1]][ip == "Sig70"]

exo <- do.call(rbind,exo)
                   
files <- files[sapply(exo[,(edsn)],function(x)grep(x,files))]
files <- file.path(exo_dir,files)

reads <- mclapply(files,create_reads,mc.cores = 6)

pbc <- sapply(reads,PBC)
depth <- sapply(reads,nreads)
ssd <- sapply(reads,SSD)


FSR <- function(reads)
{
  summ <- summary(reads)
  rF <- summ[,sum(readsF)]
  rR <- summ[,sum(readsR)]

  return(rF / (rF + rR))
}

fsr <- sapply(reads,FSR)

exo[,nreads := depth]
exo[,pbc := pbc]
exo[,ssd := ssd]
exo[,fsr := fsr]

ecoli.size <- data.table(V1 = "U00096",V2 = 4639221)

scc <- mclapply(reads,strand_cross_corr,shift = 1:300,chrom.sizes = ecoli.size,mc.cores = 6)

nsc <- sapply(scc,function(x)x[ ,max(cross.corr) / min(cross.corr)])


rl1 <- sapply(reads,function(x)readsF(x)[[1]][,mean(end - start + 1)])
rl2 <- sapply(reads,function(x)readsR(x)[[1]][,mean(end - start + 1)])
rl <- floor(.5 * (rl1 + rl2))


rsc <- mapply(function(x,rl){
  M <- x[,max(cross.corr)]
  m <- x[,min(cross.corr)]
  fr <- x[shift == rl, (cross.corr)]
  return( (M - fr) / (m - fr))},scc,rl)


exo[, readLength := rl]
exo[ , nsc := nsc]
exo[,rsc := rsc]

scc <- mapply(function(x,y)x[,sample := y],scc,exo[,(edsn)],SIMPLIFY = FALSE)
SCC <- do.call(rbind,scc)

exo[,which.max := SCC[,which.max(cross.corr),by = sample][,(V1)]]


pdf(file = "figs/for_paper/EColi_strand_cross_corr.pdf",width = 9 , height = 5)
ggplot(SCC,aes(shift,cross.corr,colour = sample))+geom_line(size = .8)+
  scale_color_brewer(palette = "Dark2")+theme_bw()+theme(legend.position = "top")
dev.off()




