
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

stats <- build_stats(regs[[2]],gr[[2]],mc = 20)

library(hexbin)
library(scales)
r <- viridis::viridis(100,option = "D")

figs_dir <- "figs/for_paper"

pdf(file = file.path(figs_dir,"ARC_vs_URCR_example_npos.pdf"),width = 5 , height = 5 )
p <- ggplot(stats , aes( ave_reads,cover_rate))+stat_binhex(bins = 50)+
  scale_fill_gradientn(colours = r,trans = "log10",
    labels = trans_format("log10",math_format(10^.x)))+xlim(0,4)+ylim(0,1)+
  theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust = 0))+
  xlab("Average read coverage (ARC)")+
  ylab("Unique read coverage rate (URCR)")
print(p)
print( p  %+% stats[ npos > 10])  ## npos > 10 
print( p  %+% stats[ npos > 30]) ## npos > 30
print( p  %+% stats[ npos > 50]) ## npos > 100
dev.off()


## here is how to do the animation


## beamer part http://tex.stackexchange.com/questions/34921/how-to-overlap-images-in-a-beamer-slide


make_plots <- function(arc,urcr,p,regs,reads,data)
{
  gr <- mapply(rbind,readsF(reads),readsR(reads),SIMPLIFY = FALSE)
  gr <- dt2gr(do.call(rbind,gr))

  data <- copy(data[ between(ave_reads,.9*arc , 1.1*arc)])
  data <- data[between(cover_rate,.9*urcr , 1.1*urcr)]

  point <- data[,dist(c(arc,urcr) - c(ave_reads,cover_rate)),by = match]
  point <- point[which.min(V1),(match)]

  setkey(data,match)
  if(length(point) > 1){
    point <- point[1]
  }
  point <- data[point]

  reg <- dt2gr(point[,2:4,with = FALSE])
  my_gr <- subsetByOverlaps(gr,reg)
 
  fwd <- my_gr[as.character(strand(my_gr)) == "+"]
  bwd <- my_gr[as.character(strand(my_gr)) == "-"]

  fwd <- ranges(fwd)
  end(fwd) <- start(fwd)
  fwd <- coverage(fwd)

  bwd <- ranges(bwd)
  start(bwd) <- end(bwd)
  bwd <- coverage(bwd)

  cover <- function(x,reg,nm){
    out <- data.table(coord = start(reg):end(reg))
    z1 <- cumsum(runLength(x))
    z2 <- c(0,runValue(x))
    out[ , tags := stepfun(z1,z2)(coord)]
    return(out)
  }

  fwd <- cover(fwd,reg)
  bwd <- cover(bwd,reg)

  M <- 1.2 * max(fwd[,max(tags)],bwd[,max(tags)])
  
  shift <- 1:200
  loc_cc <- local_strand_cross_corr(reads,reg,shift)

  plots <- list()
  bwd[,tags := -tags]
  
  dt1 <- rbind(fwd[,strand:="F"],bwd[,strand:="R"])
  plots[[1]] <- ggplot(dt1,aes(coord,tags,colour = strand))+geom_step()+
    theme_bw()+theme(legend.position = "top",plot.title = element_text(hjust  = 0))+ylim(- M , M)+
    scale_color_brewer(palette = "Set1")+xlab("Genomic position")+ylab("ChIP read counts")
 
  plots[[2]] <- p +
    geom_point(data = point,aes(x = ave_reads,y = cover_rate),fill = "orange",colour = "orange", size = 4,shape =16)

  plots[[3]] <- ggplot(loc_cc,aes(shift,cross.corr))+
    geom_point(size = 2,shape = 1)+
    geom_line(linetype = 2)+
    theme_bw()+theme(plot.title = element_text(hjust = 0))+ylab("local Strand Cross Correlation")+
    geom_smooth(method = "loess",se = FALSE,size = 2)+
    geom_abline(slope = 0,intercept = 0,linetype = 3,size = 1,colour = "grey")

  return(plots)
  
}

arcs <- c(.2,.2,.5,1,1,2,3,3.5)
urcrs <- c(.75,.5,.4,.5,.5,.25,.25,.2)

plots <- mcmapply(make_plots,arcs,urcrs,
     MoreArgs = list(p,regs[[2]],reads[[2]],stats[npos > 50]),
     SIMPLIFY = FALSE,mc.cores = 10,mc.preschedule = FALSE)

enrich <- lapply(plots,function(x)x[[2]])
cover <- lapply(plots,function(x)x[[1]])
local_scc <- lapply(plots,function(x)x[[3]])

pdf(file = "figs/for_paper/enrichment_example.pdf",width = 5,height = 5)
a <- lapply(enrich,print)
dev.off()

pdf(file = "figs/for_paper/coverage_example.pdf",width = 6,height = 4)
a <- lapply(cover,print)
dev.off()

pdf(file = "figs/for_paper/local_scc_example.pdf",width = 6,height = 4)
a <- lapply(local_scc,print)
dev.off()


