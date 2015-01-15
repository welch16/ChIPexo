

rm(list = ls())
library(devtools)
library(scales)
library(RColorBrewer)
library(data.table)

codeDr = "/p/keles/ChIPexo/volume4/ChIPexoQual"
load_all(codeDr)

mc = 8

filecodename = "ER-Rep1"
filecodename = "H3k27ac"



dataDr = "/p/keles/ChIPexo/volume3/Analysis/Carroll/human"

dataDr = "/p/keles/ChIPexo/volume3/Analysis/Ren"

load(file.path(dataDr,filecodename,"data",paste0(filecodename,"_boundRegionsTable.RData"))) 
load(file.path(dataDr,filecodename,"data",paste0(filecodename,"_depth.RData")))
load(file.path(dataDr,filecodename,"data",paste0(filecodename,"_reads_by_region.RData")))
load(file.path(dataDr,filecodename,"data",paste0(filecodename,"_regions.RData")))

system.time(
fwd_npos <- mclapply(reads_table,function(x){
  x[strand == "+",length(unique(start)),by=region]},mc.cores = mc)
)

system.time(
bwd_npos <- mclapply(reads_table,function(x){
  x[strand == "-",length(unique(end)),by=region]},mc.cores = mc)
)

regions_length = mclapply(regions,length,mc.cores = mc)
depths = mcmapply(depth_from_reads,reads_table,regions_length,mc.cores = mc)


merge_npos <- function(f_pos,r_pos,reg_length)
{
  npos = data.table(region = 1:reg_length,
    f_npos = as.integer(0),
    r_npos = as.integer(0))
  setkey(npos,region)
  npos[f_pos[,.(region)],f_npos:= f_pos[,.(V1)]]
  npos[r_pos[,.(region)],r_npos:= r_pos[,.(V1)]]
  npos[,n_pos:= f_npos + r_npos]
  return(npos)
}

npos_list = mcmapply(merge_npos,fwd_npos,bwd_npos,regions_length,
  SIMPLIFY=FALSE,mc.cores = mc)

f_npos = mclapply(npos_list,function(x)x[["f_npos"]],mc.cores = mc)
r_npos = mclapply(npos_list,function(x)x[["r_npos"]],mc.cores = mc)
n_pos = mclapply(npos_list,function(x)x[["n_pos"]],mc.cores = mc)

chr_depths = depths

widths= mclapply(regions,width,mc.cores= mc)


lowerBounds = c(1:5,10,15,20,25,30,35,40,45,50,75,100,125,150,200,250,500,750)

filter_npos = lapply(lowerBounds,function(x){
    mcmapply(.fn_filter,n_pos,chr_depths,MoreArgs = list(lower = x),
     SIMPLIFY = FALSE,mc.cores = mc)})

depth_width = mcmapply("/",depths,widths,SIMPLIFY=FALSE,mc.cores = mc)


filter_dep_width = lapply(lowerBounds,function(x){
  mcmapply(.fn_filter,depth_width,chr_depths,MoreArgs = list(lower = x),SIMPLIFY=FALSE,mc.cores = mc)})


tab = mcmapply(function(x,y,z){
    data.table(bound = x,npos = do.call(c,y),dw = do.call(c,z))
             },lowerBounds,filter_npos,filter_dep_width,MoreArgs = list(),SIMPLIFY = FALSE,mc.cores = mc)

tab = do.call(rbind,tab)
tab$bound = factor(tab$bound)

rf = colorRampPalette(rev(brewer.pal(11,"Spectral")))
r = rf(16)

x11()

p = ggplot(tab,aes(npos,dw))+stat_binhex(bins = 70)+
    facet_wrap(~bound,ncol =4)+
    scale_fill_gradientn(colours =r,trans='log10')+
    theme(legend.position = "top")+
    scale_y_continuous(limits = c(0,10))
p


TT = tab[,mean(dw),by=list(npos,bound)]
TT$npos = factor(TT$npos)

ggplot(TT[as.numeric(as.character(npos))<11],aes(bound,V1))+geom_path(aes(group=npos,colour=npos))+scale_y_log10()

x11()
filter_regions_plot(c(1:5,10,15,20),chr_measurement = depth_width,
  chr_depths = n_pos,mc=8,measurement_label = "",log=TRUE)
