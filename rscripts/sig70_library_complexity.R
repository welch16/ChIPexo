p
rm(list = ls())

library(GenomicAlignments)
library(data.table)
library(ggplot2)
library(ChIPUtils)

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

save(exo,file = "data/for_paper/sig70_summary.RData")


pdf(file = "figs/for_paper/EColi_strand_cross_corr.pdf",width = 9 , height = 5)
ggplot(SCC,aes(shift,cross.corr,colour = sample))+geom_line(size = .8)+
  scale_color_brewer(palette = "Dark2")+theme_bw()+theme(legend.position = "top")
dev.off()


create_regions <- function(reads,lower,rangesOnly=TRUE)
{
  stopifnot(class(reads) == "GRanges")
  stopifnot(lower > 0)
  cover <- coverage(reads)
  islands <- slice(cover,lower = lower,rangesOnly = rangesOnly)
  islands <- as(islands,"GRanges")
  return(islands)
}

gr <- lapply(reads,function(x)dt2gr(rbind(readsF(x)[[1]],readsR(x)[[1]])))
regs <- lapply(gr,create_regions,lower  = 1)

build_stats <- function(region,reads)
{
  ## fix formats and stuff
  ov <- findOverlaps(region,reads)
  reads <- gr2dt(reads)
  w <- width(region)    
  region <- gr2dt(region)
  region[ , width := w]
  region[, match := paste0(seqnames,":",start,"_",end)]
  reads[  subjectHits(ov), match := region[queryHits(ov), (match)] ]  
  reads[,strand := ifelse(strand == "+", "F","R")]

  ## get base statistics
  f <- reads[,sum(strand == "F"),by = match]
  setnames(f,names(f),c("match","f"))
  setkey(f,match)
  r <- reads[,sum(strand == "R"),by = match]
  setkey(r,match)
  setnames(r,names(r),c("match","r"))
  f_uniq <- reads[strand == "F",length(unique(start)),by = match]
  setnames(f_uniq,names(f_uniq),c("match","f_pos"))
  setkey(f_uniq,match)
  r_uniq <- reads[strand == "R",length(unique(end)),by = match]
  setnames(r_uniq,names(r_uniq),c("match","r_pos"))
  setkey(r_uniq,match)

  ## merge statistics
  stats <- merge(region,f,by = "match",allow.cartesian = TRUE)
  stats <- merge(stats,r,by = "match",allow.cartesian = TRUE)
  stats <- merge(stats,f_uniq,by = "match",allow.cartesian = TRUE,all = TRUE)
  stats <- merge(stats,r_uniq,by = "match",allow.cartesian = TRUE, all = TRUE)
  stats[is.na(f_pos), f_pos := 0]
  stats[is.na(r_pos), r_pos := 0]

  ## calculate composite stats
  stats[ , depth := f + r]
  stats[ , npos := f_pos + r_pos]
  stats[ , ave_reads := depth / width]
  stats[ , cover_rate := npos / depth]
  stats[ , fsr := f / (f + r)]

  stats[ , M := as.numeric(NA)]
  stats[ , A := as.numeric(NA)]

  stats[f > 0 & r > 0, M := log2(f * r) - 2 * log2(width)]
  stats[f > 0 & r > 0, A := log2( f/ r)]

  stats[ , strand := NULL]

  return(stats)
}

stats <- mcmapply(build_stats,regs,gr,SIMPLIFY = FALSE,mc.cores = 6)

aux <- mapply(function(x,y)x[,edsn := y],stats,exo[,(edsn)],SIMPLIFY =FALSE)
aux <- do.call(rbind,aux)
aux[ , edsn := factor(edsn , levels  = c(1311,1317, 931 , 1314,1320 , 933))]
aux[ , edsn := plyr::mapvalues(edsn , from = c(1311,1317, 931 , 1314,1320 , 933),
          to = c("Rif0_rep1","Rif20_rep1","Aerobic_rep1","Rif0_rep2","Rif20_rep2","Aerobic_rep2"))]
         

pdf(file = "figs/for_paper/Sig70_number_unique_positions_perSample.pdf")
ggplot(aux[npos > 10] , aes(npos))+  geom_histogram()+scale_x_log10()+facet_wrap( ~  edsn)
dev.off()

library(hexbin)
library(scales)
r <- viridis::viridis(100,option = "D")

pdf(file = "figs/for_paper/Sig70_enrichment.pdf",width = 9 , height = 7 )
p <- ggplot(aux , aes( ave_reads,cover_rate))+stat_binhex(bins = 50)+
  facet_wrap( ~ edsn)+scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+xlim(0,7)+ylim(0,.6)+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Average read coverage")+
  ylab("Unique read coverage rate")
print(p + ggtitle("A"))
print(p + ggtitle("All regions"))
print( p  %+% aux[ npos > 10] + ggtitle("Npos > 10")) 
print( p  %+% aux[ npos > 30] + ggtitle("Npos > 30")) 
dev.off()

enrichment <- list( data = aux , plot = p)
save(enrichment , file = "data/for_paper/enrich_plot.RData")


## signal-to-noise

### filter to regions where we can calcualte the local scc

pdf(file = "figs/for_paper/Sig70_MA_plots.pdf",width = 9 , height = 7 )
p <-ggplot(aux[f > 0 & r > 0 ] , aes( M , A))+stat_binhex(bins = 50)+
  facet_wrap( ~ edsn)+scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))
print(p + ggtitle("All regions"))
print( p  %+% aux[ npos > 10] + ggtitle("Npos > 10")) 
print( p  %+% aux[ npos > 30] + ggtitle("Npos > 30")) 
dev.off()


filter_stats <- lapply(stats,function(x)x[ f > 0 & r > 0])

local_scc <- function(stat , reads , min_npos,max_npos, nsample, shift = 1:300)
{

  stat <- copy(stat[ between(npos, min_npos,max_npos)])
  message(nrow(stat))
  if(nrow(stat) > nsample){
    stat <- stat[sample(match , nsample)]
  }

  regions <- dt2gr(stat[,2:4,with = FALSE])
  regions <- split(regions,start(regions))
  names(regions) <- stat[,(match)]

  regions <- as.list(regions)
  
  scc <- mclapply(regions,function(x){
    out <- local_strand_cross_corr(reads,x,shift)
    return(out)},mc.cores = 20)

  return(scc)
}


strata_high <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 100, max_npos = Inf,nsample = 100 ,shift = 1:150),SIMPLIFY = FALSE)

strata_med1 <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 75, max_npos = 100,nsample = 100 ,shift = 1:150),SIMPLIFY = FALSE)

strata_med2 <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 50, max_npos = 75,nsample = 100 ,shift = 1:150),SIMPLIFY = FALSE)

strata_low <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 30, max_npos = 50,nsample = 100 ,shift = 1:150),SIMPLIFY = FALSE)

strata_low2 <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 15, max_npos = 30,nsample = 100 ,shift = 1:150),SIMPLIFY = FALSE)

strata_low3 <- mapply(local_scc,filter_stats,reads,
      MoreArgs = list(min_npos = 10, max_npos = 15,nsample = 100 ,shift = 1:150),SIMPLIFY = FALSE)


local_scc_data <- list(strata_high, strata_med1,strata_med2,strata_low,strata_low2,strata_low3)
save(local_scc_data , file = "data/for_paper/sig70_local_scc_by_strata.RData")

## get it as a big table
strati <- c("high","med1","med2","low1","low2","low3")
load(file = "data/for_paper/sig70_local_scc_by_strata.RData")

names(local_scc_data) <- strati
local_scc_data <- lapply(local_scc_data,function(x,nm){
  names(x) <- nm
  return(x)},exo[,(edsn)])

local_scc_data <- lapply(local_scc_data,
  function(x){
    out <- lapply(x,function(y){
    scc <- mapply(function(curve,name)curve[,match := name],y,names(y),SIMPLIFY = FALSE)
    scc <- do.call(rbind,scc)
    return(scc)
    })
    out <- mapply(function(scc,name)scc[,edsn := name],out,names(out),SIMPLIFY = FALSE)
    out <- do.call(rbind,out)
    return(out)
  })

local_scc_data <- mapply(function(x,y)x[,strata := y],local_scc_data,strati,SIMPLIFY = FALSE)
local_scc_data <- do.call(rbind,local_scc_data)

noise <- function(shift,cross.corr)
{
  if(all( is.na(cross.corr))){
    out <- Inf
  }else{
    mod <- loess(cross.corr ~ shift)
    out <- mod$s
  }
  return(out)
}

cc_max <- function(shift,cross.corr)
{
  if(all(is.na(cross.corr))){
    out <- Inf
  }else{
    
    mod <- loess(cross.corr ~ shift)
    out <- max(predict(mod))
  }
  return(out)
}
 

ss <- local_scc_data[,noise(shift,cross.corr), by = .(match, edsn,strata)]
nsc <- local_scc_data[,cc_max(shift,cross.corr), by = .(match, edsn, strata)]
setnames(ss , names(ss), c("match","edsn","strata","noise"))
setnames(nsc , names(nsc), c("match","edsn","strata","max"))

nsc[ , edsn := NULL]
nsc[ , strata := NULL]

dat <- merge(ss,nsc,by = "match",allow.cartesian = TRUE)

dat[ , nsc := max / noise ]

dat[ , strata := factor(strata, levels = rev(strati) ) ]
dat[ , strata := plyr::mapvalues(strata ,
         from = c(
           "high",
           "med1",
           "med2",
           "low1",
           "low2",
           "low3"),
         to = c(
           "(100,Inf)",
           "(75,100)",
           "(50,75)",
           "(30,50)",
           "(15,30)",
           "(10,15)"             
           ))]


strata_plots <- list()
strata_plots[[1]] <- ggplot(dat , aes( strata , noise,colour = edsn))+geom_boxplot()+facet_grid( . ~ edsn) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local SCC fit noise")+ylim(0,.25)
strata_plots[[2]] <- ggplot(dat , aes( strata , max,colour = edsn))+geom_boxplot()+facet_grid( . ~ edsn) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Loess max")+geom_abline( slope = 0 , intercept = 0, linetype = 2)+
  ylim(0,.5)
strata_plots[[3]] <- ggplot(dat , aes( strata , nsc,colour = edsn))+geom_boxplot()+facet_grid( . ~ edsn) +
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local NSC")+ylim(-1,6)+geom_abline( slope = 0 , intercept = 0, linetype = 2)+ggtitle("B")


save(strata_plots,file = "data/for_paper/Local_SCC_indicator_by_strata.RData")

pdf(file = "figs/for_paper/Local_SCC_indicator_by_strata.pdf",width = 9,height = 4)
u <- lapply(strata_plots,print)
dev.off()

library(grid)
library(gridExtra)

pdf(file = "figs/for_paper/Local_SCC_all.pdf",width = 12,height = 14)
grid.arrange(strata_plots[[1]]+ggtitle("A"),
             strata_plots[[2]]+ggtitle("B"),
             strata_plots[[3]]+ggtitle("C"),nrow = 3)
dev.off()

plots <- list()
plots[[1]] <- ggplot(dat , aes( edsn , noise , colour = edsn))+geom_boxplot()+
  facet_grid( . ~ edsn ,drop = TRUE , space = "free",scales = "free_x")+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0) ,
        axis.ticks.x = element_blank(),axis.text.x = element_text(size = 0))+
  xlab("")+ylab("Local SCC fit noise")+ylim(0,.18)
plots[[2]] <- ggplot(dat , aes( edsn , max , colour = edsn))+geom_boxplot()+
  facet_grid( . ~ edsn ,drop = TRUE , space = "free",scales = "free_x")+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0) ,
        axis.ticks.x = element_blank(),axis.text.x = element_text(size = 0))+
  xlab("")+ylab("Loess max")+ylim(-.05, .4)+geom_abline(slope = 0,intercept = 0,linetype = 2)
plots[[3]] <- ggplot(dat , aes( edsn , nsc , colour = edsn))+geom_boxplot()+
  facet_grid( . ~ edsn ,drop = TRUE , space = "free",scales = "free_x")+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", plot.title = element_text(hjust = 0) ,
        axis.ticks.x = element_blank(),axis.text.x = element_text(size = 0))+
  xlab("")+ylab("Local NSC")+ylim(-1, 3)+geom_abline(slope = 0,intercept = 0,linetype = 2)

save(plots, file = "data/for_paper/Local_SCC_indicators.RData")

pdf(file = "figs/for_paper/Local_SCC_indicators.pdf",width = 5 , height = 6)
u <- lapply(plots,print)
dev.off()


pdf(file = "figs/for_paper/Local_SCC_indicators_all.pdf",width = 5,height = 14)
grid.arrange(plots[[1]]+ggtitle("A"),
             plots[[2]]+ggtitle("B"),
             plots[[3]]+ggtitle("C"),nrow = 3)
dev.off()

