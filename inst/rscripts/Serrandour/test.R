
## Comments:

#  1- can we check if the quality score improves or not the reads

rm(list = ls())
library(devtools)
library(scales)
library(RColorBrewer)
library(data.table)


codeDr = "/p/keles/ChIPexo/volume4/ChIPexoQual"

dataDr = "/p/keles/ChIPexo/volume3/SerrandourData/human"

## dataDr = "/p/keles/ChIPexo/volume3/RenData/BAMfiles"


load_all(codeDr)

##dataDr = "/p/keles/ChIPexo/volume3/LandickData/ChIPexo"
mc =  24
exofiles = paste0("ERR3369",c(33,50,58),".bam")

exofile = file.path(dataDr,"edsn1312_042814_qc.sorted.bam" )
exofile = file.path(dataDr,exofiles[1])

system.time(a <- bound_analysis(exofile,mc))






# Cell line: MCF7
# TF : ER1

#param = NULL
#param = ScanBamParam( what = "mapq")

chipseqfiles = paste0("ERR3369",c(41,36,52,54,59,48,53),".bam")
names(exofiles) = paste0("rep",c("A","B","C"))
names(chipseqfiles) = c(paste0("rep", rep(c("A","B","C"),each=2),"_",rep(1:2,3)),"Input")

exofiles = paste0("AY",552:554,".R1.sort.bam")


reads = readGAlignmentsFromBam(file.path(dataDr,exofiles[1]),param = param)
depth = length(reads)

## Ren H3k27ac depth 29,599,796
## Ren H3k4me3  28,794,319
##       rep2 31,818,368

## Serrandour repA 9,289,835
##            repB 11,041,833 

## Landick data edsn1310 3,909,669


gr = as(reads,"GRanges")
grF = subset(gr,strand == "+")
grR = subset(gr,strand == "-")

# experiment to define the optimal region
lowerBounds = c(1,2,3,4,5,10,15,20,25,30,35,40,45,50,75,100,125,150,200,250,500,750)
regions = create_regions(gr,1)

#regions = mclapply(lowerBounds,function(x)create_regions(gr,x),mc.cores = mc)

# the idea is to study with basic properties of the islands which is the best
# strategy to build them

## nRegions = data.table(bounds = lowerBounds,
##   nr = sapply(regions,function(x)sum(sapply(x,length))))

grF = sort(grF)
grR = sort(grR)
grF = split(grF,seqnames(grF))
grR = split(grR,seqnames(grR))


all_chr = names(seqlengths(gr))
fwd_overlaps = mcmapply(reads_overlaps,all_chr,regions,grF,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
bwd_overlaps = mcmapply(reads_overlaps,all_chr,regions,grR,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

fwd_reads = mcmapply(extract_reads,all_chr,grF,fwd_overlaps,regions, MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
bwd_reads = mcmapply(extract_reads,all_chr,grR,bwd_overlaps,regions, MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)
subset_reads = mcmapply(function(x,y)rbind(x,y),fwd_reads,bwd_reads,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

chr_Nregions = mclapply(regions,length,mc.cores = mc)
fwd_depths = mcmapply(depth_from_reads,fwd_reads,chr_Nregions,SIMPLIFY=FALSE,mc.cores = mc)
bwd_depths = mcmapply(depth_from_reads,bwd_reads,chr_Nregions,SIMPLIFY=FALSE,mc.cores = mc)
chr_depths = mcmapply("+",fwd_depths,bwd_depths,SIMPLIFY=FALSE,mc.cores= mc)

fwd_strand_ratio = mcmapply("/",fwd_depths,chr_depths,SIMPLIFY=FALSE,mc.cores = mc)
labels = mcmapply(function(x,y)ifelse(x >0,ifelse(y > 0,"both","fwd"),"bwd"),fwd_depths,bwd_depths,SIMPLIFY = FALSE,mc.cores = mc)
depth_width_ratio = mcmapply(function(chr_regions,chr_depth){
  return(chr_depth/width(chr_regions))}
  ,regions,chr_depths,MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)

M_values = mcmapply(function(chr_regions,fwd_depth,bwd_depth){
  return(fwd_depth * bwd_depth / width(chr_regions)^2)},regions,fwd_depths,bwd_depths,
  MoreArgs = list(),SIMPLIFY=FALSE,mc.cores = mc)
A_values = mcmapply("/",fwd_depths,bwd_depths,MoreArgs= list(),
  SIMPLIFY=FALSE,mc.cores =mc)

M_values = mcmapply(log10,M_values,SIMPLIFY=FALSE,mc.cores = mc)
A_values = mcmapply(log10,A_values,SIMPLIFY=FALSE,mc.cores=mc)


plots = list()
plots[[1]] = filter_regions_plot(lowerBounds,chr_depths,fwd_strand_ratio,"fwd/(fwd + bwd)")
plots[[2]] = filter_regions_plot(lowerBounds,chr_depths,depth_width_ratio,"depth/width",TRUE)
plots[[3]] = filter_label_plot(lowerBounds,chr_depths,labels)
plots[[4]] = filter_MA_plot(lowerBounds,chr_depths,M_values,A_values)

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

filtered_depth_width_ratio =list( do.call(c,as.list(filtered_depth_width_ratio)))

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
    scale_x_continuous(limits = c(0,200))

histograms[["width"]]=  ggplot(regions_table,aes(width))+geom_histogram(aes(y = ..density..),binwidth=30)+
    scale_x_continuous(limits = c(0,1000))

x11()

ggplot(regions_table[label == "both"],aes(prob))+geom_density()+geom_vline(xintercept=0.5,linetype=2)

ggplot(regions_table,aes(depth_width_ratio))+geom_histogram(aes(y=..density..))+
    scale_x_continuous(limits = c(0,3))

ggplot(regions_table[label == "both"],aes(diff))+geom_density()+scale_x_continuous(limits = c(-250,250))


histograms[[1]]+geom_vline(xintercept =  quantile(regions_table$nr_pos,probs = c(0,.25,.5,.75,.95,1)),colour = I("red"))

summary(regions_table$width)

library(hexbin)
library(classInt)

names(regions_table)


intervals =  classIntervals(regions_table$nr_pos,style = "fixed",fixedBreaks = quantile(regions_table$nr_pos,probs = c(0,.05,.25,.5,.75,.95,1)))

regions_table[,nr_pos_lab := cut(regions_table$nr_pos,intervals$brks)]

ggplot(regions_table[f > 0 & r > 0 & width > 0 , ],aes(1/sqrt(2)*log2(f*r/width^2),log2(f/r)))+geom_point()+geom_smooth(aes(group = nr_pos_lab,colour = nr_pos_lab),se=FALSE,na.rm =TRUE,method ="loess",size=1)+scale_color_brewer(palette = "Set1")+theme(legend.position = "bottom")+facet_grid(.~nr_pos_lab)


x11()
ggplot(regions_table,aes(nr_pos,diff))+geom_point()


x11()

rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r = rf(16)


ggplot(regions_table[f > 0 & r > 0 & width > 0 , ],aes(log2(f/r),diff))+stat_binhex(nrbins=80)+geom_smooth(aes(group = nr_pos_lab,colour = nr_pos_lab),se=FALSE,na.rm =TRUE,method ="loess",size=1)+scale_color_brewer(palette = "Set1")+theme(legend.position = "bottom")+facet_grid(.~nr_pos_lab)+ scale_fill_gradientn(colours = r,trans='log10')


x11()

ggplot(regions_table[f > 0 & r > 0 & width > 0 , ],aes(diff,colour = nr_pos_lab))+geom_density()+scale_color_brewer(palette = "Dark2")+theme(legend.position = "bottom")+facet_grid(.~nr_pos_lab)+scale_x_continuous(limits = 300*c(-1,1))+geom_vline(xintercept = 0,linetype =2,size=1)


x11()
ggplot(regions_table[f > 0 & r > 0 & width > 0 , ],aes(log2(f/r),diff))+stat_binhex(nrbins=60)+geom_smooth(se=FALSE,na.rm =TRUE,method ="loess",size=1)+scale_color_brewer(palette = "Set1")+theme(legend.position = "bottom")+ scale_fill_gradientn(colours = r,trans='log10')






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
   scale_y_discrete(breaks = c(seq(0,maxdepth,by=25),maxdepth) )+
   theme(legend.position = "top",axis.text.x = element_text(angle = 90))+xlab("number of positions per island")+
    ylab("number of reads per island")+ggtitle(set)
  return(p)
}

x11()
positions_reads_map(exofiles[1],regions_table[nr_pos_lab != "(184,794]"])

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


