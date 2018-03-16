rm(list = ls())

library(ChIPexoQual)
library(magrittr)

dr = "/p/keles/ChIPexo/volume4"

files = list.files(dr,recursive = TRUE,full.names = TRUE)
files = files[grep("carroll",files)]
files = files[grep("bam",files)]
files = files[grep("mouse",files)]
files = files[grep("sort",files)]
files = files[grep("bai",files,invert = TRUE)]
files = files[grep("txt",files,invert = TRUE)]

files = files[grep("chipseq",files,invert = TRUE)]

library(parallel)

options("mc.cores" = 22)

reads = files %>% mclapply(readGAlignments,param = NULL)
reads = reads %>% mclapply(as,"GRanges")
names(reads) = paste0("Rep",c(3,1,2))

reads = reads[paste0("Rep",seq_along(reads))]


## peakfiles = list.files(dr,full.names = TRUE,recursive = TRUE,pattern = "peaks")
## peakfiles = peakfiles[grep("mouse",peakfiles)]

## library(readr)
## library(dplyr)
## library(data.table)

## peaks = mclapply(peakfiles,read_delim,delim = " ",col_names = FALSE)xt
## peaks = lapply(peaks,function(x){
##     x = x %>% select(X1,X2,X3)
##     setnames(x,names(x),c("seqnames","start","end"))
##     x = as.data.table(x)
##   return(ChIPUtils::dt2gr(x))})
## names(peaks) = paste0("Rep",c(3,1,2))

fimofiles = list.files(dr,recursive = TRUE,full.names =TRUE)
fimofiles = fimofiles[grep("fimo.txt",fimofiles)]
fimofiles = fimofiles[grep("carroll",fimofiles)]
fimofiles = fimofiles[grepl("FoxA1",fimofiles) | grepl("FOXA1",fimofiles)]

fimofiles = fimofiles[c(1:9)]

fimo = lapply(fimofiles,readr::read_tsv)

names(fimo) = paste0("Rep",seq_along(fimo))

maxQval = 1e-2
window_length = 20
sm = 1

## reverse_cover = reads %>% mclapply(function(x){
##     subset(x,as.character(strand(x)) == "-") %>% coverage } )

## forward_cover = reads %>% mclapply(function(x){
##     subset(x,as.character(strand(x)) == "+") %>% coverage } )

profile_table <- function(fimo,fwd,bwd,depth,repl,maxQval,windowLength,smooth)
{


    regions = fimo %>% dplyr::filter(`q-value` <= maxQval) %>%
        tidyr::separate(`sequence name`,into  = c("chr","rStart"),sep ="\\:") %>%
        tidyr::separate(rStart,into = c("rstart","rend"),sep = "\\-")%>%
        mutate(rstart =as.numeric(rstart) , rend = as.numeric(rend)) 

    
    ## regions = fimo_repl %>% filter(`q-value` <= maxQval) %>%
    ##     tidyr::separate(`sequence name`,into = c("chr","rStart"),sep = "\\:") %>%
    ##     tidyr::separate(rStart , into = c("rstart","rend"),sep = "-") %>%
    ##     mutate(rstart =as.numeric(rstart) , rend = as.numeric(rend)) 

    fwdranges = regions %>% dplyr::filter(strand == "+") %>%
        mutate(start = start + rstart, stop  = start ) %>% DataFrame %>%
        makeGRangesFromDataFrame
    revranges = regions %>% dplyr::filter(strand == "-") %>%
        mutate(start = stop + rstart, stop  = start ) %>% DataFrame %>%
        makeGRangesFromDataFrame

    fwdranges = resize(fwdranges,2*window_length + 1,fix = "center")
    revranges = resize(revranges,2*window_length + 1,fix = "center")
    
    coord = seq( - window_length, window_length)
    fwdcounts = fwd[fwdranges] %>% lapply(as.vector)
    bwdcounts = bwd[revranges] %>% lapply(as.vector)

    fwdcounts = fwdcounts %>% lapply(function(x)tibble(coord,counts = x)) %>%
        bind_rows %>% mutate(strand = "+")
    bwdcounts = bwdcounts %>% lapply(function(x)tibble(coord,counts = x)) %>%
        bind_rows %>% mutate(strand = "-")

    rbind(fwdcounts,bwdcounts) %>% mutate(replicate = repl,signal = 1e9 * counts / depth)
}


profile_table2 <- function(fimo,fwd,bwd,depth,repl,maxPval,windowLength,smooth)
{


    regions = fimo %>% dplyr::filter(`p-value` <= maxPval) %>%
        tidyr::separate(`sequence name`,into  = c("chr","rStart"),sep ="\\:") %>%
        tidyr::separate(rStart,into = c("rstart","rend"),sep = "\\-")%>%
        mutate(rstart =as.numeric(rstart) , rend = as.numeric(rend)) 

    
    ## regions = fimo_repl %>% filter(`q-value` <= maxQval) %>%
    ##     tidyr::separate(`sequence name`,into = c("chr","rStart"),sep = "\\:") %>%
    ##     tidyr::separate(rStart , into = c("rstart","rend"),sep = "-") %>%
    ##     mutate(rstart =as.numeric(rstart) , rend = as.numeric(rend)) 

    fwdranges = regions %>% dplyr::filter(strand == "+") %>%
        mutate(start = start + rstart, stop  = start,strand = "*" ) %>% DataFrame %>%
        makeGRangesFromDataFrame
    revranges = regions %>% dplyr::filter(strand == "-") %>%
        mutate(start = stop + rstart, stop  = start,strand = "*" ) %>% DataFrame %>%
        makeGRangesFromDataFrame

    fwdranges = resize(fwdranges,2*window_length + 1,fix = "center")
    revranges = resize(revranges,2*window_length + 1,fix = "center")
    
    coord = seq( - window_length, window_length)
    fwdcounts = fwd[fwdranges] %>% lapply(as.vector)
    bwdcounts = bwd[revranges] %>% lapply(as.vector)

    fwdcounts = fwdcounts %>% lapply(function(x)tibble(coord,counts = x)) %>%
        bind_rows %>% mutate(strand = "+")
    bwdcounts = bwdcounts %>% lapply(function(x)tibble(coord,counts = x)) %>%
        bind_rows %>% mutate(strand = "-")

    rbind(fwdcounts,bwdcounts) %>% mutate(replicate = repl,signal = 1e9 * counts / depth)
}


FFF = fimo %>% map(function(x){x %>% dplyr::filter(`p-value` <= 1e-4) %>%
        tidyr::separate(`sequence name`,into  = c("chr","rStart"),sep ="\\:") %>%
        tidyr::separate(rStart,into = c("rstart","rend"),sep = "\\-")%>%
        mutate(rstart =as.numeric(rstart) , rend = as.numeric(rend)) %>%
        dplyr::filter(strand == "+") %>%
        mutate(start = start + rstart, stop  = start,strand = "*" ) %>% DataFrame %>%
        makeGRangesFromDataFrame
})


nreads <- function(fimoR , freads, rreads)
{
  fimoR = fimoR %>% resize(2 * window_length + 1)
  out = tibble(f = countOverlaps(fimoR,freads),r = countOverlaps(fimoR,rreads))
  out %>% colSums
}    

mapply(nreads,FFF[7:9],ff,rr)

mapply(nreads,FFF[c(2,3,1)],ff,rr)

mapply(nreads,FFF[c(5,7,4)],ff,rr)



fwd = reads %>% lapply(S4Vectors::subset,strand == "+")
bwd = reads %>% lapply(S4Vectors::subset,strand == "-")

fwd_cover = fwd %>% lapply(resize,width = 1) %>% lapply(coverage)
bwd_cover = bwd %>% lapply(resize,width = 1) %>% lapply(coverage)

depth = reads %>% lapply(length)


maxQval = 5e-2
window_length = 20
sm = 1


library(tidyverse)

library(ggplot2)


fimo %>% map(function(x)x %>% filter(`p-value` <= 0.05))



profiles1F = mcmapply(profile_table2,fimo[7:9],fwd_cover,fwd_cover,depth,c("Rep1","Rep2","Rep3"),
                      maxQval,window_length,sm ,SIMPLIFY = FALSE) %>% bind_rows


profiles1R = mcmapply(profile_table2,fimo[7:9],bwd_cover,bwd_cover,depth,c("Rep1","Rep2","Rep3"),
                      maxQval,window_length,sm ,SIMPLIFY = FALSE) %>% bind_rows




write_tsv(bind_rows(profiles1F %>% dplyr::mutate(Set = "Forward"),
                    profiles1R %>% dplyr::mutate(Set = "Reverse")),
          path = "data/figures/fig4/fig4_FoxA1_profiles_full_peaks.tsv")


profiles2F = mcmapply(profile_table2,fimo[c(2,3,1)],fwd_cover,fwd_cover,depth,c("Rep1","Rep2","Rep3"),
       maxQval,window_length,sm ,SIMPLIFY = FALSE) %>% bind_rows 

profiles2R = mcmapply(profile_table2,fimo[c(2,3,1)],bwd_cover,bwd_cover,depth,c("Rep1","Rep2","Rep3"),
       maxQval,window_length,sm ,SIMPLIFY = FALSE) %>% bind_rows 


write_tsv(bind_rows(profiles2F %>% dplyr::mutate(Set = "Forward"),
                    profiles2R %>% dplyr::mutate(Set = "Reverse")),
          path = "data/figures/fig4/fig4_FoxA1_profiles_full.tsv")



write_tsv(bind_rows(profiles2F,profiles2R),path = "data/figures/fig4/fig4_FoxA1_profiles_full.tsv")


profiles2 = mcmapply(profile_table2,fimo[c(2,3,1)],fwd_cover,fwd_cover,depth,c("Rep1","Rep2","Rep3"),
       5e-2,window_length,sm ,SIMPLIFY = FALSE) %>% bind_rows
profiles22 = mcmapply(profile_table2,fimo[c(2,3,1)],bwd_cover,bwd_cover,depth,c("Rep1","Rep2","Rep3"),
       5e-2,window_length,sm ,SIMPLIFY = FALSE) %>% bind_rows


## profiles3 = mcmapply(profile_table,fimo[c(5,6,4)],fwd_cover,bwd_cover,depth,c("Rep1","Rep2","Rep3"),
##        maxQval,window_length,sm ,SIMPLIFY = FALSE) %>% bind_rows





p1 = profiles1 %>% mutate(strand = ifelse(strand == "+","F","R")) %>% group_by(coord,replicate,strand) %>%
    summarize(signal = mean(counts)) %>%
    ggplot(aes(coord,signal,colour = strand))+geom_line()+
    theme(legend.position = "top")+facet_grid( replicate  ~ .,scales = "free")+ggtitle("All peaks")+
    scale_color_brewer(palette = "Set1")+scale_x_continuous(breaks = seq(-window_length,window_length,by = 10))


p2 = profiles2 %>% mutate(strand = ifelse(strand == "+","F","R")) %>% group_by(coord,replicate,strand) %>%
    summarize(signal = mean(counts)) %>%
    ggplot(aes(coord,signal,colour = strand))+geom_line()+
    theme(legend.position = "top")+facet_grid( replicate  ~ .,scales = "free")+ggtitle("ChIP islands")+
    scale_color_brewer(palette = "Set1")+scale_x_continuous(breaks = seq(-window_length,window_length,by = 10))


p22 = profiles22 %>% mutate(strand = ifelse(strand == "+","F","R")) %>% group_by(coord,replicate,strand) %>%
    summarize(signal = mean(counts)) %>%
    ggplot(aes(coord,signal,colour = strand))+geom_line()+
    theme(legend.position = "top")+facet_grid( replicate  ~ .,scales = "free")+ggtitle("ChIP islands")+
    scale_color_brewer(palette = "Set1")+scale_x_continuous(breaks = seq(-window_length,window_length,by = 10))


p3 = profiles3 %>% mutate(strand = ifelse(strand == "+","F","R")) %>% group_by(coord,replicate,strand) %>%
    summarize(signal = mean(counts)) %>%
    ggplot(aes(coord,signal,colour = strand))+geom_line()+
    theme(legend.position = "top")+facet_grid( replicate  ~ .,scales = "free")+ggtitle("ChIP islands")+
    scale_color_brewer(palette = "Set1")+scale_x_continuous(breaks = seq(-window_length,window_length,by = 10))




## > 
## > prop = bind_rows( profiles2 %>% mutate(set = "F"), profiles22 %>% mutate(set = "R"))
## > prop %>% filter(strand == "+")
## # A tibble: 576,706 × 6
##    coord counts strand replicate   signal   set
##    <int>  <int>  <chr>     <chr>    <dbl> <chr>
## 1    -20      1      +      Rep1 45.02383     F
## 2    -19      0      +      Rep1  0.00000     F
## 3    -18      0      +      Rep1  0.00000     F
## 4    -17      0      +      Rep1  0.00000     F
## 5    -16      0      +      Rep1  0.00000     F
## 6    -15      0      +      Rep1  0.00000     F
## 7    -14      0      +      Rep1  0.00000     F
## 8    -13      0      +      Rep1  0.00000     F
## 9    -12      1      +      Rep1 45.02383     F
## 10   -11      1      +      Rep1 45.02383     F
## # ... with 576,696 more rows
## > prop %>% filter(strand == "+") %>% group_by(coord,set , replicate) %>% summarize(signal = mean(counts))
## Source: local data frame [246 x 4]
## Groups: coord, set [?]

##    coord   set replicate    signal
##    <int> <chr>     <chr>     <dbl>
## 1    -20     F      Rep1 0.6892430
## 2    -20     F      Rep2 0.4689076
## 3    -20     F      Rep3 1.0347170
## 4    -20     R      Rep1 0.3419655
## 5    -20     R      Rep2 0.2672269
## 6    -20     R      Rep3 0.4483019
## 7    -19     F      Rep1 0.6936698
## 8    -19     F      Rep2 0.4621849
## 9    -19     F      Rep3 0.7818868
## 10   -19     R      Rep1 0.3519256
## # ... with 236 more rows
## > prop %>% filter(strand == "+") %>% group_by(coord,set , replicate) %>% summarize(signal = mean(counts)) %>%
## + ggplot(aes(coord,signal,colour = set)) +geom_line()
## > ggplot(aes(coord,signal,colour = set)) +geom_line()+facet_grid( replicate ~ ., scale = "free")
## Error: ggplot2 doesn't know how to deal with data of class uneval
## > ggplot(aes(coord,signal,colour = set)) +geom_line()+facet_grid( replicate ~ ., scales = "free")
## Error: ggplot2 doesn't know how to deal with data of class uneval
## > prop %>% filter(strand == "+") %>% group_by(coord,set , replicate) %>% summarize(signal = mean(counts)) %>% 
## + ggplot(aes(coord,signal,colour = set)) +geom_line()+facet_grid( replicate ~ ., scaleS = "free")
## Error in facet_grid(replicate ~ ., scaleS = "free") : 
##   unused argument (scaleS = "free")
## > prop %>% filter(strand == "+") %>% group_by(coord,set , replicate) %>% summarize(signal = mean(counts)) %>% 
## + ggplot(aes(coord,signal,colour = set)) +geom_line()+facet_grid( replicate ~ ., scales = "free")
## > p = prop %>% filter(strand == "+") %>% group_by(coord,set , replicate) %>% summarize(signal = mean(counts)) %>% 
## + ggplot(aes(coord,signal,colour = set)) +geom_line()+facet_grid( replicate ~ ., scales = "free")
## > p + theme_bw()+scale_color_brewer(palette = "Set1")
## > ggsave("article_copy.png",width = 7,height = 5)
## > p2 = p + theme_bw()+scale_color_brewer(palette = "Set1",name = "")
## > ggsave(p2,filename = "article_copy.png",width = 7,height = 5)
## > 
