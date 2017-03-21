#!/usr/bin/env Rscript

library(optparse,quietly = TRUE)

optList = list(
    make_option("--outdr", action = "store_true",type = "character",
                default = "figs/figuresV2/fig5",
                help = "Directory where all the figures are stored."),
    make_option("--fig.width",action = "store_true",type = "numeric",default = 7,
                help = "Width of the figure composed with all panels"),
    make_option("--line.width",action = "store_true",type = "numeric",default = .8,
                help = "Plot line.width")    
)

opt = parse_args(OptionParser(option_list = optList))

library(grid)
library(gridExtra)
library(tidyverse)
library(magrittr)
library(viridis)
library(scales)

indr = "data/figures/fig5"

files = list.files(indr,full.names = TRUE,pattern = "tsv")

files = files[grep("scores",files)]

files1 = files[grep("M.tsv",files,invert = TRUE)]
files2 = files[grep("M.tsv",files,invert = FALSE)]

scores = lapply(files1,read_tsv)

scores = mapply(function(x,y)x %>% mutate(samp = y),scores,basename(files1),SIMPLIFY = FALSE) %>%
    bind_rows
