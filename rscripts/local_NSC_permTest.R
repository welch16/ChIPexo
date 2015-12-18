
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

mc <- 16

dr <- "/p/keles/ChIPexo/volume3/CarrollData/mouse"
files <- list.files(dr)

files <- files[grep("bam",files)]
files <- files[grep("bai",files,invert = TRUE)]

files <- file.path(dr,files)

rep1 <- files[2]
rep2 <- files[3]
rep3 <- files[1]

rep1_reads <- create_reads(rep1)

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

gr <- do.call(rbind,mapply(rbind,readsF(rep1_reads),readsR(rep1_reads),
  SIMPLIFY = FALSE))                           

regions <- create_regions(dt2gr(gr),lower = 1)

regions <- dt2gr(regions)




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

stats <- build_stats(regions,dt2gr(gr),mc = mc)
## > nrow(stats)
## [1] 8472248

## > lapply(stats,summary)
## $width
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   34.00   36.00   36.00   42.88   36.00 3390.00 

## $f
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##     0.00     0.00     1.00     1.31     1.00 34610.00 

## $r
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##     0.00     0.00     1.00     1.31     1.00 54250.00 

## $f_pos
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##    0.0000    0.0000    1.0000    0.8314    1.0000 1701.0000 

## $r_pos
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
##    0.0000    0.0000    1.0000    0.8314    1.0000 1128.0000 

## $depth
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##     1.00     1.00     1.00     2.62     2.00 88860.00 

## $npos
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##    1.000    1.000    1.000    1.663    1.000 2829.000 

## $ave_reads
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##   0.0263   0.0278   0.0278   0.0466   0.0556 439.9000 

## $cover_rate
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.002138 0.500000 1.000000 0.810900 1.000000 1.000000 

## $fsr
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.0000  0.5000  0.5001  1.0000  1.0000 

## $M
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##     -13     -11     -10     -10      -9      15 7456785 

## $A
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
##      -9      -1       0       0       1       8 7456785 


## only keep regions where local-SCC can be calculated
stats <- stats[ f > 0 & r > 0]
stats <- stats[seqnames != "chrM"]

## > nrow(stats)
## [1] 1015439

## > lapply(stats,summary)
## $width
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   35.00   49.00   62.00   77.42   78.00 3390.00 

## $f
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##     1.00     1.00     2.00     4.84     3.00 34610.00 

## $r
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##     1.00     1.00     2.00     4.86     3.00 54250.00 

## $f_pos
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##    1.000    1.000    1.000    2.757    2.000 1701.000 

## $r_pos
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##    1.000    1.000    1.000    2.754    2.000 1128.000 

## $depth
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     2.0     2.0     3.0     9.7     6.0 88860.0 

## $npos
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##    2.000    2.000    2.000    5.511    3.000 2829.000 

## $ave_reads
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##   0.0274   0.0441   0.0588   0.0793   0.0833 439.9000 

## $cover_rate
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.002138 0.600000 0.684200 0.750600 1.000000 1.000000 

## $fsr
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.001613 0.333300 0.500000 0.500000 0.666700 0.995300 

## $M
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -12.720 -11.130 -10.420 -10.150  -9.462  15.490 

## $A
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## -9.27400 -1.00000  0.00000 -0.00027  1.00000  7.73500 

setkey(stats,match)

get_candidate_region <- function(mat,stats)dt2gr(stats[mat,2:4,with = FALSE])

cover_plot <- function(reg,reads,perm = FALSE)
{

  gr <- do.call(rbind,mapply(rbind,readsF(reads),readsR(reads),SIMPLIFY =FALSE))
  my_gr <- subsetByOverlaps(dt2gr(gr),reg)
  my_gr <- resize(my_gr,1)

  if(perm)strand(my_gr) <- sample(as.character(strand(my_gr)))

  my_gr <- split(my_gr,as.character(strand(my_gr)))
  fwd <- ranges(my_gr[["+"]])
  bwd <- ranges(my_gr[["-"]])

  fwd <- coverage(fwd)
  bwd <- coverage(bwd)

  fwdDT <- data.table(coord = start(reg):end(reg),tags = 0,strand = "F")
  bwdDT <- data.table(coord = start(reg):end(reg),tags = 0,strand = "R")

  fill_cover <- function(DT,cover)
  {
    z1 <- cumsum(runLength(cover))
    z2 <- c(0,runValue(cover))
    DT[,tags := stepfun(z1,z2)(coord)]
    return(DT)
  }

  fwdDT <- fill_cover(fwdDT,fwd)
  bwdDT <- fill_cover(bwdDT,bwd)

  bwdDT[,tags := -tags]

  DT <- rbind(fwdDT,bwdDT)

  out <- ggplot(DT,aes(coord,tags,colour = strand))+
    geom_step()+theme_bw()+
    theme(legend.position = "top",plot.title = element_text(hjust  = 0))+
      scale_color_brewer(palette = "Set1")+
    xlab("Genomic position")+ylab("ChIP read counts")

  return(out)

}

stats <- stats[width > 150]
## > nrow(stats)
## [1] 64165

## > lapply(stats,summary)

## $width
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     151     179     227     269     317    3390 

## $f
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##     1.00    10.00    19.00    43.46    43.00 34610.00 

## $r
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##     1.00    10.00    19.00    43.89    44.00 54250.00 

## $f_pos
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##     1.0     7.0    12.0    21.6    26.0  1701.0 

## $r_pos
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    1.00    7.00   12.00   21.55   26.00 1128.00 

## $depth
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##     5.00    21.00    39.00    87.35    85.00 88860.00 

## $npos
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##    5.00   14.00   24.00   43.16   51.00 2829.00 

## $ave_reads
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
##   0.0305   0.1079   0.1683   0.2598   0.2848 439.9000 

## $cover_rate
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.002138 0.548400 0.625000 0.626600 0.701300 1.000000 

## $fsr
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 0.001613 0.413800 0.500000 0.500200 0.586200 0.995300 

## $M
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## -12.720  -8.808  -7.384  -7.169  -5.763  15.490 

## $A
##      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
## -9.274000 -0.502500  0.000000 -0.001624  0.502500  7.735000 



M <- stats[npos >= 500,(match)][1:3]


pdf(height = 8,width = 6)
for(m1 in M){

## between(npos,30,50)][,(match)][7]

reg1 <- get_candidate_region(m1,stats)


truth <- local_strand_cross_corr(rep1_reads,reg1,shift = 1:150,perm = FALSE)

N <- 1000
system.time(
perms <- mclapply(1:N,function(i)
                  local_strand_cross_corr(rep1_reads,reg1,shift = 1:150,perm = TRUE),
                  mc.cores = mc,mc.preschedule = TRUE)
)

perms <- mapply(function(x,y)x[,perm := paste0("M",y)],perms,1:N,SIMPLIFY = FALSE)
perms <- do.call(rbind,perms)
perms[,perm := factor(perm)]


## medians <- perms[,median(cross.corr),by = shift]
## setnames(medians,names(medians),c("shift","cross.corr"))

means <- perms[,mean(cross.corr),by = shift]
setnames(means,names(means),c("shift","cross.corr"))


p <- cover_plot(reg1,rep1_reads)


truth[,summary(loess(cross.corr ~ shift))]
means[,summary(loess(cross.corr ~ shift))]


grid.arrange(p+ggtitle(m1),
ggplot(truth,aes(shift,cross.corr))+geom_point(shape = 1)+
  geom_smooth(method = "loess",se = FALSE)+
  geom_point(data = means,colour = "red",shape = 1)+
  geom_line(size = .1,linetype = 2)+
  geom_line(colour = "red",data = means,size = .1,linetype = 2)+
  geom_smooth(method = "loess",colour =  "red",data = means,se = FALSE)+
  geom_abline(slope = 0,intercept = 0,linetype = 2)+
  theme_bw()+geom_vline(xintercept = 35,colour = "blue",linetype = 1,size = .1),
             nrow = 2)
}
dev.off()

