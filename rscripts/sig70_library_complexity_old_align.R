
rm(list = ls())

library(GenomicAlignments)
library(data.table)
library(ggplot2)
library(ChIPUtils)

exo_dir <- "/p/keles/ChIPexo/volume6/sigma70_old_alignment"

files <- list.files(exo_dir)
files <- files[grep("R1",files)]
files <- files[grep("sort",files)]
files <- file.path(exo_dir,files)
files <- files[grep("bai",files,invert = TRUE)]

reads <- create_reads(files)
summary(reads)

pbc <- PBC(reads)
depth <- nreads(reads)
ssd <- SSD(reads)


FSR <- function(reads)
{
  summ <- summary(reads)
  rF <- summ[,sum(readsF)]
  rR <- summ[,sum(readsR)]

  return(rF / (rF + rR))
}

fsr <- FSR(reads)

exo <- data.table("sig70")
exo[,nreads := depth]
exo[,pbc := pbc]
exo[,ssd := ssd]
exo[,fsr := fsr]

ecoli.size <- data.table(V1 = names(readsF(reads)),V2 = 4639221)



scc <- strand_cross_corr(reads,shift = 1:300,chrom.sizes = ecoli.size)
nsc <- scc[ ,max(cross.corr) / min(cross.corr)]

rl1 <- readsF(reads)[[1]][,mean(end - start + 1)]
rl2 <- readsR(reads)[[1]][,mean(end - start + 1)]
rl <- floor(.5 * (rl1 + rl2))

exo[, readLength := rl]
exo[ , nsc := nsc]


create_regions <- function(reads,lower,rangesOnly=TRUE)
{
  stopifnot(class(reads) == "GRanges")
  stopifnot(lower > 0)
  cover <- coverage(reads)
  islands <- slice(cover,lower = lower,rangesOnly = rangesOnly)
  islands <- as(islands,"GRanges")
  return(islands)
}

gr <- dt2gr(rbind(readsF(reads)[[1]],readsR(reads)[[1]]))
regs <- create_regions(gr,lower  = 1)

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

stats <-build_stats(regs,gr)


library(reshape2)


forward_strand_ratio_plot <- function(stat,probs = c(0,.25,.5,.75,1),values = 1:750, mc = 8)
{
  stopifnot(is.numeric(values))
  stopifnot(length(probs) <= 9)
  stopifnot(length(values) > 1)
  to_use <- stat[,.(depth, prob)]
  rows <- mclapply(values,function(i)to_use[depth > i,quantile(prob, probs)],mc.cores = mc)
  rows <- mapply(function(x,i)data.table(depth = i,quantiles = probs,fsr = x ),rows,values,SIMPLIFY = FALSE)
  rows <- do.call(rbind,rows)
  rows[,quantiles := factor(quantiles,levels = sort(probs))]
  p <- ggplot(rows,aes(depth , fsr , colour = quantiles))+geom_line(size = 1)+
    scale_color_brewer(name = "quantiles",palette = "Set1")+ylim(0,1)+xlim(1,max(values))+theme_bw()+
    theme(legend.position = "top")+ylab("fwd strand ratio")+xlab("least amount of fragments in region")
  return(p)
}

label_plot <- function(stat,values = 1:750,mc = 8,prop = FALSE)
{
  stopifnot(is.numeric(values))
  stopifnot(length(values) > 1)

  to_use <- stat[,.(depth,label)]
  setkey(to_use,label)
  fwd <- to_use["fwd"]
  fwd <- do.call(c,mclapply(values,function(i)nrow(fwd[depth > i]),mc.cores = mc))
  bwd <- to_use["bwd"]
  bwd <- do.call(c,mclapply(values,function(i)nrow(bwd[depth > i]),mc.cores = mc))
  both <- to_use["both"]
  both <- do.call(c,mclapply(values,function(i)nrow(both[depth > i]),mc.cores = mc))
  dt <- data.table(fwd,bwd,both)  
  if(prop){
    rs <- rowSums(dt)
    dt <- dt / rs
  }
  ord <- c("both","fwd","bwd")
  setcolorder(dt,ord)
  dt <- cbind(values,dt)
  dt <- melt(dt,id.vars = "values")
  dt[ ,variable := factor(variable, levels = ord)]
  p <- ggplot(dt , aes(values,value, fill = variable))+geom_bar(stat="identity")+
    scale_fill_brewer(name = "Strand composition",palette = "Pastel1")+xlim(1,max(values))+theme_bw()+
    theme(legend.position = "top")+
    xlab("least amount of fragments in region")
  if(prop){
    p <- p + ylab("Cummulative proportion")+ylim(0,1)
  }else{
    p <- p + ylab("Cummulative nr. regions")+scale_y_log10(label = trans_format('log10',math_format(10^.x)))
  }

  return(p)
}


stats[, prob := fsr]
stats[,label := ifelse( f > 0 & r > 0  , "both",ifelse(f > 0,"fwd","bwd"))]
pdf(file = file.path(figs_dir,"strand_imbalance_analysis.pdf"))
forward_strand_ratio_plot(stats,probs = c(.1,.25,.5,.75,.9))
label_plot(stats)
label_plot(stats,prop = TRUE)
forward_strand_ratio_plot(stats[f > 0 & r > 0],probs = c(.1,.25,.5,.75,.9))+ggtitle("Conditional on f,r >0")
dev.off()



figs_dir <- "figs/old_data_qc"

pdf(file = file.path(figs_dir,"basic_qc_plots.pdf"))
ggplot(scc,aes(shift,cross.corr))+geom_line()+ggtitle("SCC")
ggplot(stats , aes(npos))+  geom_histogram()+scale_x_log10()+ggtitle("All npos")
ggplot(stats[npos > 10] , aes(npos))+  geom_histogram()+scale_x_log10()+ggtitle("Npos > 10")
dev.off()

library(hexbin)
library(scales)
r <- viridis::viridis(100,option = "D")


pdf(file = file.path(figs_dir,"Sig70_enrichment.pdf"),width = 9 , height = 7 )
p <- ggplot(stats , aes( ave_reads,cover_rate))+stat_binhex(bins = 50)+
  scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+xlim(0,7)+ylim(0,.6)+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Average read coverage")+
  ylab("Unique read coverage rate")
print(p + ggtitle("A"))
print(p + ggtitle("All regions"))
print( p  %+% stats[ npos > 10] + ggtitle("A"))  ## npos > 10 
print( p  %+% stats[ npos > 30] + ggtitle("B")) ## npos > 30
print( p  %+% stats[ npos > 50] + ggtitle("B")) ## npos > 50
dev.off()

## signal-to-noise

### filter to regions where we can calcualte the local scc

stats <- stats[ f > 0 & r > 0]


pdf(file = file.path(figs_dir,"Sig70_MA_plots.pdf"),width = 9 , height = 7 )
p <-ggplot(stats , aes( M , A))+stat_binhex(bins = 50)+
  scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))+ylim(-10,10)
print(p + ggtitle("All regions"))
print( p  %+% stats[ npos > 10] + ggtitle("Npos > 10")) 
print( p  %+% stats[ npos > 30] + ggtitle("Npos > 30")) 
dev.off()


pdf(file = file.path(figs_dir,"npos_vs_depth_by_region.pdf"))      
p <-ggplot(stats , aes( npos , depth))+stat_binhex(bins = 100)+
  scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlim(0,100)+ylim(0,1e4)+geom_smooth(method = "loess",se =FALSE)+
  xlab("number of unique positions per region")+ylab("nr of reads per region")
print(p)
dev.off()


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


strata_high <- local_scc(stats,reads,min_npos = 150, max_npos = Inf,nsample = 100 ,shift = 1:150)
strata_med1 <- local_scc(stats,reads,min_npos = 75, max_npos = 150,nsample = 100 ,shift = 1:150)
strata_med2 <- local_scc(stats,reads,min_npos = 50, max_npos = 75,nsample = 100 ,shift = 1:150)
strata_low <- local_scc(stats,reads,min_npos = 30, max_npos = 50,nsample = 100 ,shift = 1:150)
strata_low2 <- local_scc(stats,reads,min_npos = 15, max_npos = 30,nsample = 100 ,shift = 1:150)
strata_low3 <- local_scc(stats,reads,min_npos = 10, max_npos = 15,nsample = 100 ,shift = 1:150)

local_scc_data <- list(strata_high, strata_med1,strata_med2,strata_low,strata_low2,strata_low3)


## get it as a big table
strati <- c("high","med1","med2","low1","low2","low3")

local_scc_data <- lapply(local_scc_data,function(x){
  names(x) <- paste0("M",1:length(x))
  return(x)})

local_scc_data <- lapply(local_scc_data,                         
  function(x){
    out <- mapply(function(curve,name)curve[,match := name],x,names(x),SIMPLIFY = FALSE)
    scc <- do.call(rbind,out)
    return(scc)
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
 

ss <- local_scc_data[,noise(shift,cross.corr), by = .(match,strata)]
nsc <- local_scc_data[,cc_max(shift,cross.corr), by = .(match, strata)]
setnames(ss , names(ss), c("match","strata","noise"))
setnames(nsc , names(nsc), c("match","strata","max"))

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
           "(150,Inf)",
           "(75,150)",
           "(50,75)",
           "(30,50)",
           "(15,30)",
           "(10,15)"             
           ))]



            

strata_plots <- list()
strata_plots[[1]] <- ggplot(dat , aes( strata , noise,colour = strata))+geom_boxplot()+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local SCC fit noise")+ylim(0,.2)
strata_plots[[2]] <- ggplot(dat , aes( strata , max,colour = strata))+geom_boxplot()+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Loess max")+geom_abline( slope = 0 , intercept = 0, linetype = 2)+
  ylim(0,.3)
strata_plots[[3]] <- ggplot(dat , aes( strata , nsc,colour = strata))+geom_boxplot()+
  scale_color_brewer(palette = "Dark2")+theme_bw()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90),plot.title = element_text(hjust = 0))+
  xlab("")+ylab("Local NSC")+ylim(-.5,2)+geom_abline( slope = 0 , intercept = 0, linetype = 2)+ggtitle("B")


pdf(file = file.path(figs_dir,"Local_SCC_indicator_by_strata.pdf"),width = 5,height = 4)
u <- lapply(strata_plots,print)
dev.off()





### plot regions


plot_regions <- function(stat , reads , min_npos,max_npos, nsample, shift = 1:300)
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
  
  plots <- mclapply(regions,function(x){
    out <- segvis_to_go(reads,x)
    return(out)},mc.cores = 20)

  return(plots)
}

segvis_to_go <- function(reads,region)
{
  fwd <- dt2gr(readsF(reads)[[1]])
  end(fwd) <- start(fwd)
  bwd <- dt2gr(readsR(reads)[[1]])
  start(bwd) <- end(bwd)
  
  fwd <- subsetByOverlaps(fwd ,region)
  bwd <- subsetByOverlaps(bwd ,region)

  fwd <- coverage(ranges(fwd))
  bwd <- coverage(ranges(bwd))

  dt <- data.table( coord = start(region):end(region))
  dt[ , fwd := stepfun(cumsum(runLength(fwd)),c(0,runValue(fwd)))(coord)]
  dt[ , bwd := -stepfun(cumsum(runLength(bwd)),c(0,runValue(bwd)))(coord)]

  dt <- melt(dt, id.vars = "coord")

  out <- ggplot(dt , aes(coord ,value , colour = variable ))+geom_line()+
    scale_color_brewer(palette = "Set1",name = "")+theme(legend.position = "none")+
    theme_bw()
  return(out)
  
}


strata_high <- plot_regions(stats,reads,min_npos = 150, max_npos = Inf,nsample = 100 ,shift = 1:150)
strata_med1 <- plot_regions(stats,reads,min_npos = 75, max_npos = 150,nsample = 100 ,shift = 1:150)
strata_med2 <- plot_regions(stats,reads,min_npos = 50, max_npos = 75,nsample = 100 ,shift = 1:150)
strata_low <- plot_regions(stats,reads,min_npos = 30, max_npos = 50,nsample = 100 ,shift = 1:150)
strata_low2 <- plot_regions(stats,reads,min_npos = 15, max_npos = 30,nsample = 100 ,shift = 1:150)
strata_low3 <- plot_regions(stats,reads,min_npos = 10, max_npos = 15,nsample = 100 ,shift = 1:150)


pdf(file = file.path(figs_dir,"strata_high_profiles.pdf"))
strata_high
dev.off()
    

pdf(file = file.path(figs_dir,"strata_med1_profiles.pdf"))
strata_med1
dev.off()

pdf(file = file.path(figs_dir,"strata_med2_profiles.pdf"))
strata_med2
dev.off()

pdf(file = file.path(figs_dir,"strata_low_profiles.pdf"))
strata_low
dev.off()

pdf(file = file.path(figs_dir,"strata_low2_profiles.pdf"))
strata_low2
dev.off()

pdf(file = file.path(figs_dir,"strata_low3_profiles.pdf"))
strata_low3
dev.off()




