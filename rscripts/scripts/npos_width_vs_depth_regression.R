#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    npos_width_vs_depth_regression.R - Perform a linear model with number of unique
      positions and width as covariates and depth as response.

  Arguments:


   -- infile

      A RData file with the outcome of the ChIP-exo QC pipeline run

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

infile <- args[1]
outfile <- args[2]
M <- as.numeric(args[3])
N <- as.numeric(args[4])

mc <- detectCores()

load(infile)

stats <- ext_stats[[1]]


samp <- function(stats,nreg)
{
  dt <- stats[sample(.N,nreg),.(npos,width,depth)]
  model <- lm(depth ~ 0 + npos + width , data =dt)
  out <- data.table(tidy(model))
  return(out)
}

coeff <- mclapply(1:N,function(i)samp(stats,M),mc.cores = mc)
coeff <- do.call(rbind,coeff)

save(coeff,file = outfile)
