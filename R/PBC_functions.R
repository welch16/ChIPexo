
pbc_fun <- function(tags)
{
  tab = table(tags)
  nUniq = sum(tab == 1)
  nTotal = sum(tab >= 1)
  return(round(nUniq/nTotal,4))
}

#' @title PCR bottleneck coefficient function
#' @description Calculates a vector with 3 PCR bottleneck coeff, the ones conditionated to forward (or backward) reads and the complete PBC
#' @param ranged_data, GRanges or GAlignment object
#' @param isPET, logical that asses if the data was generated from a paired end tag experiment
#' @returns A numerical vector with the 3 PCR bottleneck coefficients described above
#' @export
PBC <- function(ranged_data,isPET=FALSE)
{
  fwd_reads = subset(ranged_data,subset = strand(ranged_data) == "+")
  bwd_reads = subset(ranged_data,subset = strand(ranged_data) == "-")
  if(isPET){
    tags = list(fwd = paste0(seqnames(fwd_reads),":",start(fwd_reads),"-",end(fwd_reads)),
      bwd = paste0(seqnames(bwd_reads),":",start(bwd_reads),"-",end(bwd_reads)))
    tags[["all"]] = c(tags[["fwd"]],tags[["bwd"]])
    all_pbc = mclapply(tags,FUN = pbc_fun,mc.cores = 3)
    names(all_pbc) = c("fwd","bwd","all")    
  }else{
    tags = list(fwd = paste0(seqnames(fwd_reads),":",start(fwd_reads),"F"),
      bwd = paste0(seqnames(bwd_reads),":",end(bwd_reads),"R")
      )
    tags[["all"]] = c(tags[["fwd"]],tags[["bwd"]])
    all_pbc = mclapply(tags,FUN = pbc_fun,mc.cores = 3)
    names(all_pbc) = c("fwd","bwd","all")   
  }
  all_pbc = do.call(c,all_pbc)
  return(all_pbc)
}

pbc_to_df <- function(pbc)
{
  df = melt(pbc)
  colnames(df) = c("set","type","pbc")
  df$set = do.call(c,lapply(df$set,function(x,st)gsub("_042814","",x),df$set))
  return(df)
}
