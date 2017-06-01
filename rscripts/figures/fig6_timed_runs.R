rm(list = ls())


library(tidyverse)

indr = "data/timed_runs"
files = indr %>%
    list.files(full.names = TRUE,
               recursive = TRUE)

times = files %>%
    map( read_tsv ) %>%
    map2( files,
         .f = function(x,y)
             x %>%
             mutate(
                 samp = gsub("_time.tsv","",basename(y))
             )) %>%
    bind_rows() %>%
    filter(task != "Sampling reads") %>%    
    mutate(cores = factor(cores,levels = c(4,8,16)),
           id = as.character(seq_len(nrow(.))),
           task = factor(task,levels = unique(task) %>% rev()))

samp_dt = times %>%
    filter(cores == 4) %>%
    group_by(samp) %>%
    summarize(
        total_time = sum(time)) %>%
    arrange(total_time)

samp_lvls = samp_dt %>%
    select(samp) %>% .[[1]]

theme_set(theme_bw())

lvs = times %>% filter(cores == 4) %>%
    filter(grepl("ChIPexoQual",task)) %>%
    arrange(time) %>% select(samp) %>% .[[1]]

## times %>% 
##     ggplot(aes(samp,time,fill = task))+
##     geom_bar(stat = "identity")+
##     theme(legend.position = "top",
##           axis.text.x = element_text(hjust = .5),
##           axis.title.y = element_blank())+
##     facet_grid(cores ~ .,scales = "free_x")+coord_flip()

## times %>% filter(cores == 4) %>%
##     ggplot(aes(samp,time,fill = task))+
##     geom_bar(stat = "identity")+
##     theme(legend.position = "top",
##           axis.text.x = element_text(hjust = .5),
##           axis.title.y = element_blank())+
##     coord_flip()+ylab("Time (in seconds)")+
##     scale_fill_brewer(palette = "Set1")    

depths = tibble(
    samp = c("ER_rep2",
             "CTCF_rep2",
             "CTCF_rep1",
             "GR_IMR90",
             "TBPnexus_K562",
             "TBPexo_K562"),
    depths = c(11041833,
               20947081,
               23576694,
               47443803,
               129675001,
               114282270)) %>%
    mutate(
        depths = prettyNum(depths,big.mark = ","))

times = times %>% left_join(depths,by = "samp")



times = times %>%
    mutate(depths = ifelse(task == "Reading aligned reads",
                           "",depths))  %>%
    left_join(samp_dt,by = "samp") %>%
    mutate(
        samp = factor(samp,levels = samp_lvls))

pdf("figs/figuresV2/times_comp.pdf",height = 3,width = 5)
times %>% filter(cores == 4) %>%
    ggplot(aes(x = samp,y = time,fill = task))+
    geom_bar(stat = "identity")+
    geom_text(aes(x = samp,y = total_time + 250,label = depths),hjust = 1)+
    theme(legend.position = "top",
          axis.text.x = element_text(hjust = .5),
          axis.title.y = element_blank())+
    coord_flip()+ylab("Time (in seconds)")+
    scale_fill_brewer(palette = "Set1",name = "")+
    ylim(0,1e3)
times %>% filter(cores == 4 &
                 !grepl("Reading",task)) %>%
    mutate(samp =factor(samp,levels = lvs )) %>%
    ggplot(aes(x = samp,y = time))+
    geom_bar(stat = "identity")+
    geom_text(aes(x = samp,y = time + 145,label = depths),hjust = 1)+
    theme(legend.position = "top",
          axis.text.x = element_text(hjust = .5),
          axis.title.y = element_blank())+
    coord_flip()+ylab("Time (in seconds)")
dev.off()


## times %>% filter(cores == 4) %>%
##     ggplot(aes(x = samp,y = time,fill = task))+
##     geom_bar(stat = "identity")+
##     geom_text(aes(x = samp,y = total_time + 130,label = depths),hjust = 1)+
##     theme(legend.position = "top",
##           axis.text.x = element_text(hjust = .5),
##           axis.title.y = element_blank())+
##     coord_flip()+ylab("Time (in seconds)")+
##     scale_fill_brewer(palette = "Set1",name = "")+
##     ylim(0,1e3)
## times %>% filter(cores == 4 &
##                  !grepl("Reading",task)) %>%
##     mutate(samp =factor(samp,levels = lvs )) %>%
##     ggplot(aes(x = samp,y = time))+
##     geom_bar(stat = "identity")+
##     geom_text(aes(x = samp,y = time + 70,label = depths),hjust = 1)+
##     theme(legend.position = "top",
##           axis.text.x = element_text(hjust = .5),
##           axis.title.y = element_blank())+
##     coord_flip()+ylab("Time (in seconds)")


## > times %>% filter(grepl("ChIPexo",task)) %>% group_by(cores) %>%
## + summarize( min = min(time),
## + max = max(time))
## # A tibble: 3 x 3
##    cores     min     max
##   <fctr>   <dbl>   <dbl>
## 1      4  82.780 421.344
## 2      8 120.333 610.698
## 3     16 138.599 908.631
