
rm(list = ls())

library(data.table)
library(parallel)
library(reshape2)
library(ggplot2)
library(scales)
library(ggrepel)
library(devtools)
load_all("~/Desktop/Docs/Code/ChIPUtils")

exodir <- "data/ChIPexo_QC_runs"
exofiles <- list.files(exodir)

coeffdir <- "data/ARC_quantify"
files <- list.files(coeffdir)


load_file <- function(x){load(x);coeff}

coeffs <- lapply(file.path(coeffdir,files),load_file)

coeffs <- mapply(function(x,y){
  x$id = y
  x},coeffs,gsub(".RData","",files),SIMPLIFY = FALSE)


coeffs <- data.table(do.call(rbind,coeffs))


size_dir <- system.file("extdata","chrom.sizes",package = "ChIPUtils")
sizes <- lapply(file.path(size_dir,list.files(size_dir)),read.table)
sizes <- lapply(sizes,data.table)
names(sizes) <- gsub(".chrom.sizes","",list.files(size_dir))
sizes <- lapply(sizes,function(x)x[V1 != "chrM",sum(as.numeric(V2))])
sizes <- mapply(function(x,y){
  data.table(V1 = y, V2 = x)},sizes,names(sizes),SIMPLIFY = FALSE)
ecoli.size <- data.table(V1 = "U00096",V2 = 4639221)
sizes <- rbind(do.call(rbind,sizes),ecoli.size)
sizes <- sizes[V1 != "mm10"]


load_exo <- function(x){
  load(x)  
  return(ext_stats[["nreads"]])}

depth <- mclapply(file.path(exodir,exofiles),load_exo,mc.cores = detectCores())

depth <- data.table(id = gsub(".RData","",exofiles),
                    depth =do.call(c,depth))

coeffs <- dcast.data.table(data = coeffs[,.(id,term,estimate)],
  formula = id ~ term, value.var = "estimate" )

dt <- merge(depth,coeffs,by = "id")

dt[,organism := "Human"]
dt[grepl("mouse",id),organism := "Mouse"]
dt[grepl("landick",id),organism := "E.Coli"]
dt[grepl("S2",id),organism := "D.Melanogaster"]
dt[grepl("embryo",id),organism := "D.Melanogaster"]
dt[ , seq := "ChIP-exo"]
dt[ grepl("nexus",id), seq := "ChIP-nexus"]
dt[,depth := depth / 1e6]

setnames(sizes,names(sizes),c("genome","genome.size"))

sizes[,organism := c("D.Melanogaster","Human","Mouse","E.Coli")]

dt <- merge(dt,sizes,by = "organism")

dt[,norm_depth := depth / genome.size]

s <- 3
pdf(file = "figs/test_quantify/arc_vs_urcr_param_against_depth.pdf",width = 12,height = 6)
ggplot(dt,aes(norm_depth,b,colour = organism,fill = organism))+
  geom_point(size = s)+
  scale_color_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "top")+
  geom_text_repel(data = dt,aes(norm_depth,b,label = id,colour = organism),show_guide = FALSE,size = 3)+
  scale_x_log10()
du <- dt[id != "venters_TBP_K562_Rep3"]
ggplot(du,aes(norm_depth,k,colour = organism,fill = organism))+
  geom_point(size = s)+
  scale_color_brewer(palette = "Set1")+scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "top")+
  scale_x_log10()+
  geom_text_repel(data = du,aes(norm_depth,k,label = id,colour = organism),show_guide = FALSE,size = 3)
dev.off()

modk <- lm( k ~  depth * organism, data = dt)
summary(modk)

## Call:
## lm(formula = k ~ depth * organism, data = dt[organism != "Mouse"])

## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.046069 -0.003801 -0.000183  0.002610  0.083401 

## Coefficients:
##                       Estimate Std. Error t value Pr(>|t|)
## (Intercept)          1.831e-02  1.544e-02   1.186    0.248
## depth                2.132e-04  6.293e-04   0.339    0.738
## organismE.Coli       3.488e-02  2.193e-02   1.591    0.126
## organismHuman        4.222e-03  1.915e-02   0.221    0.827
## depth:organismE.Coli 1.129e-05  1.528e-03   0.007    0.994
## depth:organismHuman  2.175e-04  6.487e-04   0.335    0.741

## Residual standard error: 0.02411 on 22 degrees of freedom
## Multiple R-squared:  0.4166,	Adjusted R-squared:  0.2841 
## F-statistic: 3.142 on 5 and 22 DF,  p-value: 0.02733


modb <- lm( b ~ depth*organism , data = dt)
summary(modb)

## Call:
## lm(formula = b ~ depth * organism, data = dt[organism != "Mouse"])

## Residuals:
##       Min        1Q    Median        3Q       Max 
## -0.260250 -0.038366 -0.000404  0.068970  0.183123 

## Coefficients:
##                       Estimate Std. Error t value Pr(>|t|)    
## (Intercept)           0.423157   0.076311   5.545 1.42e-05 ***
## depth                -0.001636   0.003111  -0.526  0.60424    
## organismE.Coli       -0.362288   0.108413  -3.342  0.00295 ** 
## organismHuman         0.122490   0.094653   1.294  0.20906    
## depth:organismE.Coli -0.001737   0.007552  -0.230  0.82020    
## depth:organismHuman  -0.003416   0.003207  -1.065  0.29833    
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Residual standard error: 0.1192 on 22 degrees of freedom
## Multiple R-squared:  0.7852,	Adjusted R-squared:  0.7364 
## F-statistic: 16.09 on 5 and 22 DF,  p-value: 1.041e-06
