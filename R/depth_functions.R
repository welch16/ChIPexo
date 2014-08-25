
as.GRanges <- function(x)return(as(x,"GRanges"))

logcount <- function(set,binsize)
{
  bin = create.bins(binsize,seqlengths(set))
  count = countOverlaps(bin,set)
  return(log(1 + count))
}

logcount_data.frame <- function(binsize,sets,mc.cores = mc.cores)
{
  logcounts = mclapply(sets,FUN = logcount,binsize,mc.cores = mc.cores)
  dataframe = melt(logcounts)
  colnames(dataframe) = c("log_counts","sample")
  dataframe$binSize = binsize
  return(dataframe)
}

mergeDataFrames <- function(df_list)
{
  df = do.call(rbind,df_list)
  df$sample = factor(df$sample)
  df$binSize = factor(df$binSize)
  return(df)
}

logcount_boxplot <- function(dataframe)
{
  p = ggplot(dataframe,aes(binSize,log_counts,fill = sample))+geom_boxplot()+
    theme(axis.text.x = element_text(angle = 90),legend.position = "none",
        strip.text.x = element_text(size = 8,angle = 90))+          
    facet_grid(. ~sample)
  return(p)
}


assignDepth <- function(seq.list,depth)
{
  seq.list$depth = rep(-1,nrow(seq.list))
  idx = do.call(c,lapply(1:nrow(seq.list),function(i,seq.list,depth)
    grep(seq.list[i,1],names(depth)),seq.list,depth))
  seq.list$depth[idx] = depth[idx]
  return(seq.list)
}
