
rm(list = ls())

library(ChIPQC)
library(parallel)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)

dr = "otherQC_methods/ChIPQC"
files = list.files(dr,full.names = TRUE)

load(files[1])
load(files[2])
load(files[3])


