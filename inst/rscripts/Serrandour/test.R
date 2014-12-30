
## Comments:

#  1- can we check if the quality score improves or not the reads

rm(list = ls())
library(devtools)
library(scales)
library(RColorBrewer)

codeDr = "/p/keles/ChIPexo/volume4/ChIPexoQual"
dataDr = "/p/keles/ChIPexo/volume3/SerrandourData/human"

load_all(codeDr)

# Cell line: MCF7
# TF : ER1

param = NULL  # param = ScanBamParam( what = "mapq")
mc =18

exofiles = paste0("ERR3369",c(33,50,58),".bam")
chipseqfiles = paste0("ERR3369",c(41,36,52,54,59,48,53),".bam")
names(exofiles) = paste0("rep",c("A","B","C"))
names(chipseqfiles) = c(paste0("rep", rep(c("A","B","C"),each=2),"_",rep(1:2,3)),"Input")

reads = readGAlignmentsFromBam(file.path(dataDr,exofiles[1]),param = param)
depth = length(reads)

gr = as(reads,"GRanges")
grF = subset(gr,strand == "+")
grR = subset(gr,strand == "-")

# experiment to define the optimal region
lowerBounds = c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,75,100,125,150,200,250,500,750)
regions = mclapply(lowerBounds,function(x)create_regions(gr,x),mc.cores = mc)

# the idea is to study with basic properties of the islands which is the best
# strategy to build them

nRegions = data.table(bounds = lowerBounds,
  nr = sapply(regions,function(x)sum(sapply(x,length))))

grF = sort(grF)
grR = sort(grR)
grF = split(grF,seqnames(grF))
grR = split(grR,seqnames(grR))


all_chr = names(seqlengths(gr))
fwd_overlaps = mcmapply(reads_overlaps,all_chr,regions[[1]],grF,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
bwd_overlaps = mcmapply(reads_overlaps,all_chr,regions[[1]],grR,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

extract_reads <- function(chr,chr_reads,chr_overlaps,chr_islands)
{
  # need to check this function, there is an environment issue with
  # calling data.table from the inside of the function
  chr_overlaps = .data.table.Hits(chr_overlaps)
  chr_reads = .data.table.GRanges(chr_reads)
  setkey(chr_overlaps,queryHits)
  chr_reads$ID = 1:nrow(chr_reads)
  setkey(chr_reads,ID)
  chr_islands = .data.table.GRanges(GRanges(seqnames = chr,
    ranges = chr_islands,strand = "*")) # check here, possible bug
  ov_reads = chr_reads[chr_overlaps[queryHits %in% 1:nrow(chr_islands),list(subjectHits)],]
  ov_reads[,ID := NULL]
  ov_reads[,region:=chr_overlaps[,queryHits]]
  return(ov_reads)
}

fwd_reads = mcmapply(extract_reads,all_chr,grF,fwd_overlaps,regions[[1]], MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
bwd_reads = mcmapply(extract_reads,all_chr,grR,bwd_overlaps,regions[[1]], MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

subset_reads = mcmapply(function(x,y)rbind(x,y),fwd_reads,bwd_reads,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

reads2IRanges <- function(x)IRanges(start = x$start,end = x$end)

depth_from_reads <- function(chr_reads,chr_nregions= NULL)
{
  if(is.null(chr_nregions)){    
    chr_depth = table(chr_reads$region)

  }else{
    chr_depth = table(factor(
      chr_reads$region,levels = 1:chr_nregions))     
  }
  dimnames(chr_depth) = NULL
  return(chr_depth)
}

chr_Nregions = mclapply(regions[[1]],length,mc.cores = mc)

fwd_depths = mcmapply(depth_from_reads,fwd_reads,chr_Nregions,SIMPLIFY=FALSE,mc.cores = mc)
bwd_depths = mcmapply(depth_from_reads,bwd_reads,chr_Nregions,SIMPLIFY=FALSE,mc.cores = mc)
chr_depths = mcmapply(function(x,y)x+y,fwd_depths,bwd_depths,SIMPLIFY=FALSE,mc.cores= mc)

fwd_strand_ratio = mcmapply(function(x,y)x/y,fwd_depths,chr_depths,SIMPLIFY=FALSE,mc.cores = mc)
labels = mcmapply(function(x,y)ifelse(x >0,ifelse(y > 0,"both","fwd"),"bwd"),fwd_depths,bwd_depths,SIMPLIFY = FALSE,mc.cores = mc)

depth_width_ratio = mcmapply(function(chr_regions,chr_depth){
  return(chr_depth/width(chr_regions))}
  ,regions[[1]],chr_depths,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)


fn_filter <- function(chr_ratio,chr_depth,lower)return(chr_ratio[chr_depth > lower])


filter_regions_plots <- function(lowerBounds,chr_depths,chr_measurement,measurement_label = "")
{
  filter_regions = lapply(lowerBounds,function(x){
    mcmapply(fn_filter,chr_measurement,chr_depths,MoreArgs = list(lower = x),
      SIMPLIFY = FALSE,mc.cores = mc)})           
  names(filter_regions) = lowerBounds
 
  filter_regions = mcmapply(function(x,y){
    data.table(bound = x,measurement = do.call(c,y))},lowerBounds,filter_regions,MoreArgs = list(),
    SIMPLIFY = FALSE,mc.cores = mc)

  filter_regions = do.call(rbind,filter_regions)
  filter_regions$bound = factor(filter_regions$bound)

  p1 = ggplot(filter_regions,aes(bound,measurement,colour=bound))+
    geom_boxplot(outlier.colour = alpha("black",1/50))+
    theme(legend.position = "none")+
    ylab(measurement_label)

  p2 = ggplot(filter_regions[,length(measurement),by=bound],aes(bound,V1,fill=bound))+
    geom_bar(stat = "identity")+
    scale_x_discrete(breaks = lowerBounds)+
    xlab("bound")+
    scale_y_log10(  labels=trans_format('log10', math_format(10^.x)))+
    ylab("Nr. islands")+
    theme(legend.position ="none")  
  return(list(p2,p1))
}


filter_log_regions_plots <- function(lowerBounds,chr_depths,chr_measurement,measurement_label = "")
{
  # Need to fix this function, I am re-using the depth of wider regions
  # when increasing the lower depth, i.e. the depths need to change not filter

  # I still think it is more noticeable in the highest values of lowerBounds
  filter_regions = lapply(lowerBounds,function(x){
    mcmapply(fn_filter,chr_measurement,chr_depths,MoreArgs = list(lower = x),
      SIMPLIFY = FALSE,mc.cores = mc)})           
  names(filter_regions) = lowerBounds
 
  filter_regions = mcmapply(function(x,y){
    data.table(bound = x,measurement = do.call(c,y))},lowerBounds,filter_regions,MoreArgs = list(),
    SIMPLIFY = FALSE,mc.cores = mc)

  filter_regions = do.call(rbind,filter_regions)
  filter_regions$bound = factor(filter_regions$bound)

  p1 = ggplot(filter_regions,aes(bound,measurement,colour=bound))+
    geom_boxplot(outlier.colour = alpha("black",1/50))+
    theme(legend.position = "none")+
    scale_y_log10(  labels=trans_format('log10', math_format(10^.x)))+
    ylab(measurement_label)

  return(p1)
}



plots = filter_regions_plots(lowerBounds,chr_depths,fwd_strand_ratio,"fwd/(fwd + bwd)")

plots[[3]] = filter_log_regions_plots(lowerBounds,chr_depths,depth_width_ratio,"depth/width")

filter_regions_plots <- function(lowerBounds,chr_depths,labels)
{
  filter_labels= lapply(lowerBounds,function(x){
    mcmapply(fn_filter,labels,chr_depths,MoreArgs = list(lower = x),
    SIMPLIFY = FALSE,mc.cores = mc)})           
  names(filter_labels) = lowerBounds

  filter_labels = mcmapply(function(x,y){
    data.table(bound = x,measurement = do.call(c,y))},lowerBounds,filter_labels,MoreArgs = list(),
    SIMPLIFY = FALSE,mc.cores = mc)

  filter_labels = do.call(rbind,filter_labels)
  filter_labels$bound = factor(filter_labels$bound)
  filter_labels$measurement = factor(filter_labels$measurement)
 
  p = ggplot(filter_labels,aes(bound,fill = measurement))+
    geom_bar(position = "fill")+
    scale_fill_brewer(palette = "Set1")+
    theme(legend.position = "bottom")+
    ylab("proportion by label")
  
  return(p) 
}

plots[[4]] = filter_regions_plots(lowerBounds,chr_depths,labels)

u =lapply(plots,function(x){x11(width = 6,height = 5);print(x)})
graphics.off()



# The idea for this plots is to use the forward strand ratio to get the best filter lower bound

# it seems that there is one bound after which the filter stabilizes
# in this case 15 or 20...

# also for the depth, that seems to happening, i.e. after certain bound the log10 number of regions seems to
# be stabilizing (for example from 15 to 50)..

# I think that perhaps we are removing tons of false regions between those bound values

# Summary stats:
# island statistics as start, end , width (implied in filtered_regions)
# strand_depth,depth_fwd_strand_ratio --> need to filter
# depth / width

# coverage statistics
# nr.positions , summit_position , plots , 

# remove islands with low depth

bound = 15
which_regions = mclapply(chr_depths,function(x,bound)which(x > bound),bound,mc.cores = mc)

#sum(sapply(which_regions,length))
# bound  nr.Regions
#  15    30,917
#  20    24,311
#  25    20,216

filtered_regions = mcmapply(function(chr_region,idxs){
  chr_region = GRanges(seqnames= "not",ranges = chr_region,strand = "*")
  return(ranges(chr_region[idxs])) 
},regions[[1]],which_regions,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)



npos <- function(reads,strand)
{
  if(strand == "+"){
    pos = nrow(reads[strand=="+",start[1],by=start])
  }else if(strand == "-"){
    pos = nrow(reads[strand=="-",end[1],by=end])
  }
  return(pos)
}

nrPos_from_reads <- function(chr,chr_reads,chr_which,strand)
{
  message("Calculating number of positions for ",chr," islands")
  pos = do.call(c,mclapply( chr_which,function(i,reads){
    npos(reads[region == i & strand == strand],strand)},
    chr_reads,mc.cores = mc))
  return(pos)
}

fwd_filter_npos = mapply(nrPos_from_reads,all_chr,subset_reads,which_regions,
  MoreArgs = list(strand = "+"),SIMPLIFY = FALSE)

bwd_filter_npos = mapply(nrPos_from_reads,all_chr,subset_reads,which_regions,
  MoreArgs = list(strand = "-"),SIMPLIFY = FALSE)

filter_npos = mcmapply(function(fwd,bwd)return(fwd+bwd),
  fwd_filter_npos,bwd_filter_npos,SIMPLIFY=FALSE,mc.cores = mc)


coverage_from_reads <- function(chr,chr_reads,chr_which,strand)
{
  message("Calculating ",ifelse(strand == "+","forward","backward"),
          " coverage for ",chr," islands")
  cover = mclapply( chr_which,function(i,reads,orientation){
    region_reads = reads[region == i & strand == orientation]
    iranges = reads2IRanges(region_reads)   
    return(coverage(iranges))
    },chr_reads,strand,mc.cores = mc)
  return(cover)
}

fwd_filtered_cover = mapply(coverage_from_reads,all_chr,subset_reads,which_regions,
  MoreArgs = list(strand = "+"),SIMPLIFY = FALSE)
bwd_filtered_cover = mapply(coverage_from_reads,all_chr,subset_reads,which_regions,
  MoreArgs = list(strand = "-"),SIMPLIFY = FALSE)


## Plotting functions

Rle2data.table <- function(rle_data)
{
  coord = cumsum(runLength(rle_data))
  counts = runValue(rle_data)
  return(data.table(coord=coord,counts = counts))
}


normalize.tagcounts <- function(counts,depth)return(counts*1e6 / depth)

plot_cover <- function(fwd_cover,bwd_cover,region_start,region_end,
                       depth,annot = NULL)
{
  
  fwd = Rle2data.table(fwd_cover)
  bwd = Rle2data.table(bwd_cover)

  fwd$strand = "fwd"
  bwd$strand = "bwd"

  M = max(c(fwd$count,bwd$count))*1.2
  
  fwd$count = normalize.tagcounts(fwd$count,depth)
  bwd$count = -normalize.tagcounts(bwd$count,depth)

  if(nrow(fwd) > 0 & nrow(bwd) > 0){
    colors = c("blue","red")
  }else{
    if(nrow(fwd) > 0){
      colors = "red"
    }else{
      colors = "blue"
    }
  }

  p = ggplot(rbind(fwd,bwd),aes(coord,count,colour = strand))+
    geom_step()+
    theme(legend.position = "none")+
    scale_colour_manual(values = colors)+
    geom_abline(slope = 0,intercept = 0,linetype = 2,colour = I("black"))+
    scale_y_continuous( limits = 1e6/depth*M * c(-1,1))+ylab("Normalized counts")

  return(p)
}


which.max.range <- function(cover)
{
  if(nrun(cover)==0){
    out = NULL
  }else{
    maxCount = max(cover)
    coords = cumsum(runLength(cover))
    tagCounts = runValue(cover)
    maxPos = which(tagCounts == maxCount)
    start = coords[maxPos]
    end = coords[maxPos + 1]
    out = matrix(c(start,end),ncol=2,byrow =FALSE)
  }
  return(out)      
}


fwd_filtered_maxRanges = lapply(fwd_filtered_cover,
  function(x,mc){
    return(mclapply(x,which.max.range,mc.cores=mc))
  },mc = mc)

bwd_filtered_maxRanges = lapply(bwd_filtered_cover,
  function(x,mc){
    return(mclapply(x,which.max.range,mc.cores=mc))
  },mc = mc)

summit_diff <- function(fwd_max_range,bwd_max_range,criteria = median)
{
  if(is.null(fwd_max_range) | is.null(bwd_max_range)){
    out = NA
  }else{
    out = criteria(bwd_max_range[,2]) - criteria(fwd_max_range[,1])    
  }
  return(out)
}

filtered_summit_diff = lapply(all_chr,function(chr,fwd_filtered_maxRanges,bwd_filtered_maxRanges,mc){
    return(mcmapply(summit_diff,
        fwd_filtered_maxRanges[[chr]],
        bwd_filtered_maxRanges[[chr]],MoreArgs = list(),mc.cores = mc))},
    fwd_filtered_maxRanges,bwd_filtered_maxRanges,mc)


# filter regions, depth and ratios


filtered_depths = mcmapply(function(chr_depth,idxs){
  return(chr_depth[idxs])},chr_depths,which_regions,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)

fwd_filtered_depths = mcmapply(function(chr_depth,idxs){
  return(chr_depth[idxs])},fwd_depths,which_regions,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)
bwd_filtered_depths = mcmapply(function(chr_depth,idxs){
  return(chr_depth[idxs])},bwd_depths,which_regions,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)

filtered_fwdstrand_ratio= mcmapply(function(chr_depth,idxs){
  return(chr_depth[idxs])},fwd_strand_ratio,which_regions,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)
filtered_labels = mcmapply(function(chr_depth,idxs){
  return(chr_depth[idxs])},labels,which_regions,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)



filtered_depth_width_ratio = mcmapply(function(chr_region,chr_depth){ 
    return(chr_depth/width(chr_region))
},filtered_regions,filtered_depths,MoreArgs = list(),mc.cores= mc)


chr_data.table <- function(chr,region)
{
  return(data.table(seqnames = chr,start =start(region),end =end(region),width = width(region)))
}

# Variables so far
# regions: filtered_regions
# number of positions per islads: filter_npos, fwd_filter_npos, bwd_filter_npos
# depth per islands:  filtered_depths, fwd_filtered_depths,bwd_filtered_depths
# ratios: filtered_fwdstrand_ratio filtered_width_depth_ratio
# maxRanges : fwd_filtered_maxRanges, bwd_filtered_maxRanges
# labels: filtered_labels
# coverages: fwd_filtered_cover, bwd_filtered_cover
# summit_diff: filtered_summit_diff

# regions: filtered_regions
filtered_regions = mcmapply(chr_data.table,all_chr,filtered_regions,MoreArgs=list(),SIMPLIFY=FALSE,mc.cores=mc)

add_variable <- function(region,var,name)
{
  region[[name]] = var
  return(region)
}

# number of positions per islads: filter_npos, fwd_filter_npos, bwd_filter_npos
filtered_regions= mcmapply(add_variable,
  filtered_regions,filter_npos,MoreArgs = list(name = "nr_pos"),SIMPLIFY=FALSE,mc.cores=mc)
filtered_regions= mcmapply(add_variable,
  filtered_regions,fwd_filter_npos,MoreArgs = list(name = "fwd_nr_pos"),SIMPLIFY=FALSE,mc.cores=mc)
filtered_regions= mcmapply(add_variable,
  filtered_regions,bwd_filter_npos,MoreArgs = list(name = "bwd_nr_pos"),SIMPLIFY=FALSE,mc.cores=mc)

# depth per islands:  filtered_depths, fwd_filtered_depths,bwd_filtered_depths
filtered_regions= mcmapply(add_variable,
  filtered_regions,filtered_depths,MoreArgs = list(name = "depth"),SIMPLIFY=FALSE,mc.cores=mc)
filtered_regions= mcmapply(add_variable,
  filtered_regions,fwd_filtered_depths,MoreArgs = list(name = "f"),SIMPLIFY=FALSE,mc.cores=mc)
filtered_regions= mcmapply(add_variable,
  filtered_regions,bwd_filtered_depths,MoreArgs = list(name = "r"),SIMPLIFY=FALSE,mc.cores=mc)

# ratios: filtered_fwdstrand_ratio filtered_width_depth_ratio
filtered_regions= mcmapply(add_variable,
  filtered_regions,filtered_fwdstrand_ratio,MoreArgs = list(name = "prob"),SIMPLIFY=FALSE,mc.cores=mc)
filtered_regions= mcmapply(add_variable,
  filtered_regions,filtered_depth_width_ratio,MoreArgs = list(name = "depth_width_ratio"),SIMPLIFY=FALSE,mc.cores=mc)

# labels: filtered_labels
# summit_diff: filtered_summit_diff
filtered_regions= mcmapply(add_variable,
  filtered_regions,filtered_labels,MoreArgs = list(name = "label"),SIMPLIFY=FALSE,mc.cores=mc)
filtered_regions= mcmapply(add_variable,
  filtered_regions,filtered_summit_diff,MoreArgs = list(name = "diff"),SIMPLIFY=FALSE,mc.cores=mc)

regions_table = do.call(rbind,filtered_regions)

# histogram / density for each variable

names(regions_table)

 ## "seqnames"          "start"             "end"              
 ## "width"             "nr_pos"            "fwd_nr_pos"       
 ## "bwd_nr_pos"        "depth"             "f"                
 ## "r"                 "prob"              "width_depth_ratio"
 ## "label"             "diff"             

histograms = list()
histograms[["nr_pos"]] = ggplot(regions_table,aes(nr_pos))+geom_histogram(aes(y = ..density..),binwidth=1)+
    scale_x_continuous(limits = c(0,500))

histograms[["nr_pos"]]=  ggplot(regions_table,aes(width))+geom_histogram(aes(y = ..density..),binwidth=30)+
    scale_x_continuous(limits = c(0,1000))

ggplot(regions_table[label == "both"],aes(prob))+geom_density()+geom_vline(xintercept=0.5,linetype=2)

ggplot(regions_table,aes(width_depth_ratio))+geom_histogram(aes(y=..density..))+
    scale_x_continuous(limits = c(0,3))

ggplot(regions_table[label == "both"],aes(diff))+geom_density()+scale_x_continuous(limits = c(-250,250))


summary(regions_table$width)


positions_reads_map <- function(set,islands,mp=Inf,md=Inf)
{
  df = islands
  setkey(df,depth,nr_pos)
  mat = df[,length(seqnames),by=list(depth,nr_pos)]
  setnames(mat,names(mat),c(names(mat)[1:2],"nr_islands"))
  rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
  r = rf(16)
  max_pos = mp
  max_depth = md
  maxdepth = max(mat[nr_pos <= max_pos & depth <= max_depth,list(depth)])
  p = ggplot(mat[nr_pos <= max_pos],aes(as.factor(nr_pos),as.factor(depth),fill = nr_islands))+
   geom_tile()+scale_fill_gradientn(name = "number of islands",colours = r,trans="log10")+
   scale_y_discrete(breaks = c(seq(0,maxdepth,by=100),maxdepth) )+
   theme(legend.position = "top",axis.text.x = element_text(angle = 90))+xlab("number of positions per island")+
    ylab("number of reads per island")+ggtitle(set)
  return(p)
}
 
positions_reads_map(exofiles[1],regions_table,mp=60,md=200)

## u = unlist(lapply(fwd_filtered_maxRanges,function(x){
##   return(lapply(x,nrow))}))
## v = unlist( lapply(bwd_filtered_maxRanges,function(x){
##   return(lapply(x,nrow))}))

## > table(u)
## u
##     1     2     3     4     5     6     7     8 
## 17933  5239  1437   382    98    25    11     1 
## > table(v)
## v
##     1     2     3     4     5     6     7     8 
## 18023  5245  1376   351    88    28     4     1 

# This is to check the number of maxima for each strand
# We can see that for most of the islands, they have a unique
# maxima. Thus, it may be worthwhile to check if this is
# useful to make a quality control

# So far, we are going to use median as a criteria to select
# the central point. This assumption may not be correct.




## chr = "chr10"
## i = 12
## x11()
## fwd_cover = fwd_filtered_cover[[chr]][[i]]
## bwd_cover = bwd_filtered_cover[[chr]][[i]]
## region_start= start(filtered_regions[[chr]])[i]
## region_end = end(filtered_regions[[chr]])[i]
## plot_cover(fwd_cover,bwd_cover,region_start,region_end,depth)


