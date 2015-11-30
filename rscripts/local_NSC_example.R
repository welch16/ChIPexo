
rm(list = ls())

library(reshape2)
library(ChIPUtils)
library(ggplot2)
library(data.table)
library(GenomicAlignments)
library(parallel)
library(viridis)
library(grid)
library(gridExtra)

dr <- "/p/keles/ChIPexo/volume3/CarrollData/mouse"
files <- list.files(dr)

files <- files[grep("bam",files)]
files <- files[grep("bai",files,invert = TRUE)]

files <- file.path(dr,files)

reads <- mclapply(files,create_reads,mc.cores = 3)
names(reads) <- files


bamfiles <-c("ERR336942.bam","ERR336956.bam","ERR336935.bam")
repl <- paste0("rep-",1:3)


create_regions <- function(reads,lower,rangesOnly=TRUE)
{
  stopifnot(class(reads) == "GRanges")
  stopifnot(lower > 0)
  sreads <- split(reads,seqnames(reads))
  sreads <- as.list(sreads)
  cover <- mcmapply(function(x,y)coverage(x)[[y]],sreads,names(sreads),mc.cores = detectCores())
  islands <- lapply(cover,slice,lower = lower, rangesOnly = rangesOnly)
  islands <- mapply(function(x,y)GRanges(seqnames = y,ranges = x),islands,names(islands),SIMPLIFY =FALSE)
  islands <- do.call(rbind,lapply(islands,gr2dt))
  return(islands)
}

gr <- lapply(reads,function(x){
  byChr <- mapply(rbind,readsF(x),readsR(x),SIMPLIFY =FALSE)
  out <- do.call(rbind,byChr)
  return(dt2gr(out))
})

regs <- lapply(gr,create_regions,lower  = 1)
regs <- lapply(regs,dt2gr)



build_stats <- function(region,reads,mc)
{
  ## fix formats and stuff
  stat_by_chr <- function(reg,read)
  {
     ov <- findOverlaps(reg,read)  
     read <- gr2dt(read)
     w <- width(reg)    
     region <- gr2dt(reg)
     region[ , width := w]
     region[, match := paste0(seqnames,":",start,"_",end)]
     read[  subjectHits(ov), match := region[queryHits(ov), (match)] ]  
     read[,strand := ifelse(strand == "+", "F","R")]
     read <- read[!is.na(match)]

     ## get base statistics
     f <- read[,sum(strand == "F"),by = match]
     setnames(f,names(f),c("match","f"))
     setkey(f,match)
     r <- read[,sum(strand == "R"),by = match]
     setkey(r,match)
     setnames(r,names(r),c("match","r"))
     f_uniq <- read[strand == "F",length(unique(start)),by = match]
     setnames(f_uniq,names(f_uniq),c("match","f_pos"))
     setkey(f_uniq,match)
     r_uniq <- read[strand == "R",length(unique(end)),by = match]
     setnames(r_uniq,names(r_uniq),c("match","r_pos"))
     setkey(r_uniq,match)

     ## merge statistics
     stats <- merge(region,f,by = "match",allow.cartesian = TRUE)
     stats <- merge(stats,r,by = "match",allow.cartesian = TRUE)
     stats <- merge(stats,f_uniq,by = "match",
                    allow.cartesian = TRUE,all = TRUE)
     stats <- merge(stats,r_uniq,by = "match",
                    allow.cartesian = TRUE, all = TRUE)
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

  reg <- split(region,as.character(seqnames(region)))
  read <- split(reads,as.character(seqnames(reads)))

  stat_chr <- mcmapply(stat_by_chr,reg,read,SIMPLIFY = FALSE,mc.cores = mc)

  stats <- do.call(rbind,stat_chr)

  return(stats)
}

stats <- mapply(build_stats,regs,gr,mc = 20 ,SIMPLIFY = FALSE)



aux <- mapply(function(x,y)x[,sample := y],stats,basename(names(stats)),SIMPLIFY =FALSE)
aux <- do.call(rbind,aux)
aux[, sample := plyr::mapvalues(sample , from = c("ERR336935.bam","ERR336942.bam","ERR336956.bam"),
    to = c("rep-3","rep-1","rep-2"))]

library(hexbin)
library(scales)
r <- viridis::viridis(100,option = "D")

figs_dir <- "figs/for_paper"

pdf(file = file.path(figs_dir,"FoxA1_enrichment.pdf"),width = 9 , height = 5 )
p <- ggplot(aux , aes( ave_reads,cover_rate))+stat_binhex(bins = 50)+
  facet_wrap( ~ sample)+scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+xlim(0,4)+ylim(0,1)+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Average read coverage (ARC)")+
  ylab("Unique read coverage rate (URCR)")
print(p + ggtitle("All regions"))
print( p  %+% aux[ npos > 10] + ggtitle("A"))  ## npos > 10 
print( p  %+% aux[ npos > 30] + ggtitle("B")) ## npos > 30
dev.off()

## here is how to do the animation


## beamer part http://tex.stackexchange.com/questions/34921/how-to-overlap-images-in-a-beamer-slide

## > du = data.table(ave_reads = 2, cover_rate = .25)
## > p %+% aux[sample == "rep-1"] + geom_point(data.table(ave_reads = 2,cover_rate = .25),aes(x = ave_reads,y = cover_rate),colour = "red", size = 2)
## Error: ggplot2 doesn't know how to deal with data of class uneval
## > p %+% aux[sample == "rep-1"] + geom_point(data = du,aes(x = ave_reads,y = cover_rate),colour = "red", size = 2)
## Warning message:
## Removed 67 rows containing missing values (stat_hexbin). 
## > dev.off()
## null device 
##           1 
## > p %+% aux[sample == "rep-1"] + geom_point(data = du,aes(x = ave_reads,y = cover_rate),colour = "red", size = 3,shape =2)
## Warning message:
## Removed 67 rows containing missing values (stat_hexbin). 
## > dev.off()
## null device 
##           1 
## > p %+% aux[sample == "rep-1"] + geom_point(data = du,aes(x = ave_reads,y = cover_rate),colour = "orange", size = 3,shape =1)
## Warning message:
## Removed 67 rows containing missing values (stat_hexbin). 
## >   C-c C-c
## > dev.off()
## null device 
##           1 
## >p %+% aux[sample == "rep-1"] + geom_point(data = du,aes(x = ave_reads,y = cover_rate),colour = "red", size   C-c C-c
## > 


all_stats <- stats

## candidates are high complexity regions
candidates <- stats[[2]][ npos > 500][seqnames != "chrM"]

candidates <- dt2gr(candidates[,2:4,with = FALSE])

ov <- lapply(regs,function(x)which(countOverlaps(x,candidates) > 0))

regs <- mapply(function(r,i)r[i],regs,ov,SIMPLIFY = FALSE)

## stats <- mapply(build_stats,regs,gr,mc = 1 ,SIMPLIFY = FALSE)

make_plots <- function(reg,regs,reads,nms)
{
  gr <- lapply(reads,function(x){
    byChr <- mapply(rbind,readsF(x),readsR(x),SIMPLIFY =FALSE)
    out <- do.call(rbind,byChr)
    return(dt2gr(out))
  })

  my_regs <- lapply(regs,function(x)x[countOverlaps(x,reg) > 0])

  my_gr <- mapply(subsetByOverlaps,gr,my_regs,SIMPLIFY = FALSE)

  fwd <- lapply(my_gr,function(x)x[as.character(strand(x)) == "+"])
  bwd <- lapply(my_gr,function(x)x[as.character(strand(x)) == "-"])

  fwd <- lapply(fwd,ranges)
  fwd <- lapply(fwd,function(x){
    end(x) <- start(x)
    return(x)})
  fwd <- lapply(fwd,coverage)

  bwd <- lapply(bwd,ranges)
  bwd <- lapply(bwd,function(x){
    start(x) <- end(x)
    return(x)})
  bwd <- lapply(bwd,coverage)

  simplified_regs <- lapply(my_regs,function(x){
    GRanges(seqnames = unique(as.character(seqnames(x))),
            ranges = IRanges(start = min(start(x)),
              end = max(end(x))))})
  
  fwd <- mapply(function(x,reg,nm){
    out <- data.table(coord = start(reg):end(reg))
    z1 <- cumsum(runLength(x))
    z2 <- c(0,runValue(x))
    out[ , tags := stepfun(z1,z2)(coord)]
    out[ ,sample := nm]
    return(out)
  },fwd,simplified_regs,nms,SIMPLIFY = FALSE)

  bwd <- mapply(function(x,reg,nm){
    out <- data.table(coord = start(reg):end(reg))
    z1 <- cumsum(runLength(x))
    z2 <- c(0,runValue(x))
    out[ , tags := stepfun(z1,z2)(coord)]
    out[ ,sample := nm]
    return(out)    
  },bwd,simplified_regs,nms,SIMPLIFY = FALSE)

  fwd <- do.call(rbind,fwd)
  bwd <- do.call(rbind,bwd)

  shift <- 1:200
  loc_cc <- mapply(local_strand_cross_corr,reads,simplified_regs,
    MoreArgs = list(shift),SIMPLIFY = FALSE)

  plots <- list()
  bwd[,tags := -tags]
  
  dt1 <- rbind(fwd[,strand:="F"],bwd[,strand:="R"])
  plots[[1]] <- ggplot(dt1,aes(coord,tags,colour = strand))+geom_step()+facet_grid( sample ~.)+
    theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust  = 0))+
    scale_color_brewer(palette = "Set1")+xlab("Genomic position")+ylab("ChIP read counts")

  loc_cc <- mapply(function(x,y )x[,sample := y],loc_cc,nms,SIMPLIFY = FALSE)
  dt2 <- do.call(rbind,loc_cc)

  plots[[2]] <- ggplot(dt2,aes(shift,cross.corr))+geom_point(size = 1,shape = 1)+
    geom_line(linetype = 2)+
    theme_bw()+theme(plot.title = element_text(hjust = 0))+ylab("local Strand Cross Correlation")+
    facet_grid(sample ~ .)

  return(plots)
  
}

nms <- c("rep-3","rep-1","rep-2")
Z <- make_plots(regs[[2]][1],regs,reads,nms)

vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

pdf(file = "figs/for_paper/local_SCC_example.pdf",width = 9,height = 5)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 5)))
dev.off()

pdf(file = "figs/for_paper/local_SCC_separated.pdf",width = 6 ,height =4)
print(Z[[1]])
print(Z[[2]] )
print(Z[[2]] +geom_smooth(method = "loess",se = FALSE ,size = 1 ))
dev.off()


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

 Z[[2]]$data[,noise(shift,cross.corr),by = sample]
 Z[[2]]$data[,cc_max(shift,cross.corr),by = sample]
 Z[[2]]$data[,cc_max(shift,cross.corr)/noise(shift,cross.corr),by = sample]
