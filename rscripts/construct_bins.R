#! /usr/bin/env Rscript
args <- commandArgs(TRUE)
if("--help" %in% args){
  cat("

  Name:

    construct_bins.R - Samples a SET file out of a PET file a write it in the same location as file

  Arguments:

    -- file

      A bam file, with the reads of a ChIP experiment. Can have Paired End Tags (PET) or
      Single End Tags (SET).

   -- outfile

      Name of the bam file with location where the sampled reads are going to be saved.

   -- N

      Number of reads to be sampled. If the reads are paired ended, we are going to sample
      2 * N reads instead of N to complete the pairs.

   -- isPET

      Logical value indicating if the reads are Paired End or Single End.

   -- seed

      Seed used to sample the reads.

   -- help

      Show the help file.

  Author:

    Rene Wech, Department of Statistics, University of Wisconsin - Madison
   
");q()}
