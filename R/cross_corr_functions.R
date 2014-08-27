

read_chipqc <- function(files,mc.cores)
{
  out = mclapply(files,function(x)ChIPQCsample(x,annotation = NULL,runCrossCor = TRUE),mc.cores = mc.cores)
  return(out)
}

read_spp <- function(files,mc.cores)
{
  out = mclapply(files,function(x)read.bam.tags(x),mc.cores =8)
  return(out)
}


cross_corr_spp <- function(seq_spp,range,bin)
{  
  out = lapply(seq_spp,function(x,range,bin)get.binding.characteristics(x,srange= range,bin = bin),range,bin)
  return(out)    
}


cross_corr_df <- function(seq_chipqc,seq_crossCorr_spp,seq_names,shift,mc.cores)
{
  n = length(seq_chipqc) 
  stopifnot( n == length(seq_crossCorr_spp))
  stopifnot( n == length(seq_names))
  df_list = mclapply(1:n,function(i,seq_chipqc,seq_crossCorr_spp,seq_names,shift){
    qc_crossCorr = seq_chipqc[[i]]@CrossCorrelation
    spp_crossCorr = seq_crossCorr_spp[[i]]$cross.correlation$y
    df = melt(data.frame(qc = qc_crossCorr,spp=spp_crossCorr))
    names(df) = c("method","crossCorr")
    df$shift = rep(shift,2)
    df$sample = seq_names[i]
    return(df)
  },seq_chipqc,seq_crossCorr_spp,seq_names,shift,mc.cores =mc.cores)
  return(df_list)
}

plot_template <- function(df_list,main,nrow)
{    
  df = do.call(rbind,df_list)
  df$sample = factor(df$sample)
  p =ggplot(df,aes(shift,crossCorr,colour = method))+geom_line()+
    facet_wrap(~sample,scales = "free",nrow = nrow)+
    theme(legend.position = "top")+ggtitle(main)
  return(p)
}
