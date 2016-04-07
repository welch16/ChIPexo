#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    arc_vs_urcr_regression.R - Performs a non-linear models between the ARC (average reads coverage) and
      the URCR (unique read coverage rate)

  Arguments:


   -- bamfile

      File with bam format that contains the fragments of the ChIP-exo experiment.

   -- outfile

      RData file used to store the run and estimated parameters.

   -- minNpos

      Min. number of unique position per region considered.

   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 3)

library(GenomicAlignments)
library(parallel)
library(devtools)
library(data.table)
library(broom)

load_all("~/Desktop/Docs/Code/ChIPexoQual")

bamfile <- args[1]
outfile <- args[2]
minNpos <- as.numeric(args[3])

mc <- detectCores()


exo <- create_exo_experiment(bamfile,parallel = TRUE,mc.cores = mc)
stats <- summary_stats(exo)

rm(exo)

stats <- stats[ f > 0 & r > 0 ]

stats <- stats[npos > minNpos]

model <- nls(cover_rate ~ b + k / ave_reads,data =stats,start = list(k =1,b=0))
coeff <- tidy(model)
coeff$n <- nrow(stats)

save(coeff,file = outfile)


## dr <- "/p/keles/ChIPexo/volume4/carroll_data/mouse/"

## files <- list.files(dr)
## files <- files[grep("sort",files)]
## files <- files[grep("bai",files,invert = TRUE)]

## bamfile <- file.path(dr,files[3])

## samp <- function(stats,M)
## {
##   idx <- sample( nrow(stats),M)
##   dt <- stats[idx,.(npos,depth)]
##   model <- lm(depth ~ 0 + npos, data =dt)
##   summ <- summary(model)
##   out <- data.table(tidy(model))
##   out[,term := NULL]
##   out[,R2 := summ$r.squared]  
##   return(out)

## }

## coeff <- mclapply(1:N,function(i)samp(stats,M),mc.cores = mc)
## coeff <- do.call(rbind,coeff)

