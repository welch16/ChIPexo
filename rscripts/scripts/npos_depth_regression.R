#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    npos_depth_regression.R - Perform a linear model between the number of unique positions and
      the number of fragments per region.

  Arguments:


   -- bamfile

      File with bam format that contains the fragments of the ChIP-exo experiment.

   -- outfile

      RData file used to store the run and estimated parameters.

   -- M

      Number of regions being drawn

   -- N

      Number of samples

   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}

stopifnot(length(args) == 4)

library(GenomicAlignments)
library(parallel)
library(devtools)
library(data.table)
library(broom)

load_all("~/Desktop/Docs/Code/ChIPexoQual")

bamfile <- args[1]
outfile <- args[2]
M <- as.numeric(args[3])
N <- as.numeric(args[4])

mc <- detectCores()

exo <- create_exo_experiment(bamfile,parallel = TRUE,mc.cores = mc)
stats <- summary_stats(exo)

rm(exo)

samp <- function(stats,M)
{
  idx <- sample( nrow(stats),M)
  dt <- stats[idx,.(npos,depth)]
  model <- lm(depth ~ 0 + npos, data =dt)
  summ <- summary(model)
  out <- data.table(tidy(model))
  out[,term := NULL]
  out[,R2 := summ$r.squared]  
  return(out)

}

coeff <- mclapply(1:N,function(i)samp(stats,M),mc.cores = mc)
coeff <- do.call(rbind,coeff)

save(coeff,file = outfile)
