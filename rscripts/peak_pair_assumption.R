

rm(list = ls())

library(data.table)
library(GenomicAlignments)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(ChIPUtils)

######################################################################################

## Initial parameters

mc <- detectCores()
figs_dir <- "figs/condition"

base_dir <- "/p/keles/ChIPexo/volume6/condition"

check_create <- function(dr)if(!dir.exists(dr))dir.create(dr)


figs_dir <- file.path(figs_dir,"all_together")
check_create(figs_dir)

######################################################################################

## get promoter locations table

prom_dir <- "/p/keles/ChIPexo/volume6/gold_standard/Landick"
pfiles <- list.files(prom_dir)
pfiles <- pfiles[grep(".bed",pfiles)]
promoters <- lapply(file.path(prom_dir,pfiles),read.table)
promoters <- lapply(promoters,data.table)

prom <- read.table(
    "/NO_BACKUP/KelesGroup_DongjunChung/ChIP-exo/Landick_ChIP-exo/database_regulonDB/PromoterSigma70Set.txt",
  skip=33, sep="\t")
loc.prom <- prom[,4]



promoters[[1]][,str := "R"]
promoters[[2]][,str := "F"]

sites <- do.call(rbind,promoters)
sites[, id := 1:nrow(sites)]
setnames(sites,names(sites),c("seqnames","start","end","strand","id"))

sites[,seqnames := "U00096"]
sites[,strand := ifelse(strand == "F","+","-")]

## zites <- data.table(seqnames = "U00096",start = loc.prom, end = loc.prom , strand = "*")

## ######################################################################################

out_dir <- "/p/keles/ChIPexo/volume6/peak_pair"
bs_dir <- file.path(out_dir,"binding")
check_create(bs_dir)

  load(file = file.path(bs_dir, "dpeakFit_exo_sig70aerobic.RData"))
  load( file = file.path(bs_dir, "dpeakFit_pet3_sig70aerobic.RData"))
  load( file = file.path(bs_dir, "dpeakFit_set_sig70aerobic.RData"))

## fit_comparison <- function(fit1 , fit2 , mm)
## {
##   mu1 <- lapply(fit1@optFit,function(x)x$mu)
##   mu2 <- lapply(fit2@optFit,function(x)x$mu)

##   G1 <- sapply(mu1,length)
##   G2 <- sapply(mu2,length)

##   print(identical(mu1,mu2))
##   print(table(G1))
##   print(table(G2))
##   print(table(G1 , G2))

##   min_dist <- function(m1,m2){

##     out <- sapply(m1 , function(x) min(abs(x - m2)))
##     out <- min(out)

##     return(out)

##   }
  
##   alls <- mapply(c , mu1 , mu2 ,SIMPLIFY = FALSE)
##   alls <- lapply(alls,sort)

##   alls <- lapply(alls,diff)

##   alls <- sapply(alls,mean)
  
##   dt <- data.table( G1 , G2  , mean_dist = alls, min = mapply(min_dist,mu1,mu2))
##   dt[,diff := G1 - G2 ]

##   r <- viridis::viridis(100, option = "D")
  
##   p1 <- ggplot(dt[,median(min),by = .(G1,G2)], aes(G1, G2 , fill = V1))+geom_tile()+
##     scale_fill_gradientn(colours = r,name = "")+ggtitle(mm)+xlab("common")+ylab("separate")+coord_flip()+
##     theme(legend.position = "top")
##   print(p1)
##   return(dt)

## }


## exo_out <- fit_comparison(exo_fits[[1]],exo_fits[[2]],"exo")
## pet_out <- fit_comparison(pet_fits1[[1]],pet_fits1[[2]],"pet")
## set_out <- fit_comparison(set_fits[[1]],set_fits[[2]],"set")
## dev.off()

site_fit_match <- function(fit,sites,mm)
{
  sites <- dt2gr(sites)
  end(sites) <- start(sites)
  strand(sites) <- "*"
  sites <- reduce(sites)

  mu <- fit@optMu

  peaks <- strsplit(names(mu),"_",fixed = TRUE)
  peaks <- GRanges(seqnames = sapply(peaks,function(x)x[1]),
                   ranges = IRanges(
                     start = as.numeric(sapply(peaks,function(x)x[2])),
                     end = as.numeric(sapply(peaks,function(x)x[3]))))

  message(mm)
  message("Nr. peaks: ",length(peaks))
  
  ov <- findOverlaps(peaks,sites)
  ov <- data.table(as.data.frame(ov))
  ov <- split(ov, ov[,(queryHits)])
  
  message("Nr. sites (low res): ",length(ov))
  
  out <- mclapply(names(ov),function(i,mu,sites,peaks){
    est <- mu[as.numeric(i)][[1]]
    names(est) <- NULL
    truth <- start(sites[ov[[i]][,(subjectHits)]])
    aux <- list()
    aux[["peak"]] <- peaks[as.numeric(i)]
    aux[["truth"]] <- truth
    aux[["est"]] <- est
    return(aux)},mu,sites,peaks,mc.cores = 16)

  return(out)

}

exo_match1 <- site_fit_match(exo_fits[[1]],sites,"exo1")
exo_match2 <- site_fit_match(exo_fits[[2]],sites,"exo2")
pet_match <- site_fit_match(pet_fits1[[1]],sites,"pet")
set_match1 <- site_fit_match(set_fits[[1]],sites,"set1")
set_match2 <- site_fit_match(set_fits[[2]],sites,"set2")


reso_table <- function(match,mm)
{
  reso <- lapply(match,function(x){
    out <- sapply(x$truth,function(z)min(abs(z - x$est)))
  return(out)})

  reso <- unlist(reso)
  out <- data.table(reso , seq = mm)

  return(out)
}

reso1 <- do.call(rbind,list(
  reso_table(exo_match1,"exo_common"),
  reso_table(exo_match2,"exo_sep"),
  reso_table(pet_match,"pet"),
  reso_table(set_match1,"set_common"),
  reso_table(set_match2,"set_sep")))



exo_match1 <- site_fit_match(exo_fits[[1]],zites,"exo1")
exo_match2 <- site_fit_match(exo_fits[[2]],zites,"exo2")
pet_match <- site_fit_match(pet_fits1[[1]],zites,"pet")
set_match1 <- site_fit_match(set_fits[[1]],zites,"set1")
set_match2 <- site_fit_match(set_fits[[2]],zites,"set2")

reso2 <- do.call(rbind,list(
  reso_table(exo_match1,"exo_common"),
  reso_table(exo_match2,"exo_sep"),
  reso_table(pet_match,"pet"),
  reso_table(set_match1,"set_common"),
  reso_table(set_match2,"set_sep")))


ggplot(reso1[reso < 100],aes(seq,reso))+geom_boxplot()+ylim(0,100)
dev.off()
