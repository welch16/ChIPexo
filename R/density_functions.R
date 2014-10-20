#' @title density.reads.per.strand.ratio
#' @description Given a predefined bin partition of the genome, this counts the number of reads with forward and backward strand in each bin and obtaines the ratio for the forward strand defined as (counts_F + 1)/(counts_F + counts_R + 2)
#' @param bins GRanges / GAlignments object
#' @param reads GRanges / GAlgnments, it need to contain the strand
#' @return density object
#' @export
#' @seealso \code{\link{create.bin}}

density.reads.per.strand.ratio <- function(bins,reads)
{  
  counts_F = countOverlaps(bins,subset(reads,subset = strand(reads) == "+"))
  counts_R = countOverlaps(bins,subset(reads,subset = strand(reads) == "-"))
  ratio = (counts_F + 1)/(counts_F + counts_R + 2) # used pseudo-counts to avoid zero denominator
  return(density(ratio))
}

add_chr <-  function(chr,bins)
{
  return(GRanges(seqnames = chr,ranges = bins,strand = "*"))
}

  

#' @title create.bins
#' @description Create a GRanges object given the length of the genome and the desired size for each bin
#' @param binSize Numeric value of the desired size for each bin
#' @param genomeLengths Numeric Length of the genome
#' @return GRanges object
#' @export
create.bins <-  function(binSize, genomeLengths)
{
  bins_list = lapply(genomeLengths,function(x,binSize)IRanges(start = seq(1,x,by=binSize),width = binSize),binSize)  
  gr_list = mapply(add_chr,names(genomeLengths),bins_list)
  bins = unlist(GRangesList(gr_list))
  seqlengths(bins) = genomeLengths
  return(bins)
}

#' @title resume.samples
#' @description Generates a string used to match ChIP-Seq PET and ChIP-Exo samples with same conditions
#' @param edsn Identifier of experiment
#' @param cult Cultivation number of the experiment
#' @param ip Inmunoprecipitate used in the experiment
#' @param phase Phase used, could be Exponential or Stationary
#' @param Growth condition of the experiment, could be Aerobic or Anaerobic
#' @param rif condition of the experiment, could be 0 min or 20 min
#' @param rep Replicate the experiment, usually 1 or 2
#' @param seq Type of sequencing used, Exo, PET or SET
#' @return Character with the logical structure used to filter the conditions table
#' @export
resume.samples <- function(edsn = NULL,cult = NULL,ip = NULL,phase = NULL,growth = NULL,
  rif = NULL,rep = NULL,seq = NULL)
{
  empty = TRUE
  st = ""
  if(!is.null(edsn)){st = paste0(st,ifelse(!empty,"&",""),"edsn=='",edsn,"'");empty = FALSE}
  if(!is.null(cult)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(cult),"is.na(cult)",paste0("cult==",cult)));empty = FALSE}
  if(!is.null(ip)){st = paste0(st,ifelse(!empty,"&",""),"ip=='",ip,"'");empty = FALSE}
  if(!is.null(phase)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(phase),"is.na(phase)",paste0("phase=='",phase,"'")))
    empty = FALSE}
  if(!is.null(growth)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(growth),"is.na(growth)",paste0("growth=='",growth,"'")))
    empty = FALSE}  
  if(!is.null(rif)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(rif),"is.na(rif)",paste0("rif=='",rif,"'")));empty = FALSE}
  if(!is.null(rep)){st = paste0(st,ifelse(!empty,"&",""),ifelse(is.na(rep),"is.na(rep)",paste0("rep==",rep)));empty = FALSE}
  if(!is.null(seq))st = paste0(st,ifelse(!empty,"&",""),"seq=='",seq,"'") 
  return(st)
}

multi.df <- function(binsizes,exo.sets,pet.sets)
{
  bins = mclapply(binsizes,FUN = create.bins ,seqlengths(exo.sets[[1]]),mc.cores =length(binsizes))

  exo.densities = suppressWarnings(mclapply(bins,function(b,seq.sets){
    z = list()
    z[[1]] = density.reads.per.strand.ratio(b,seq.sets[[1]])
    z[[2]] = density.reads.per.strand.ratio(b,seq.sets[[2]])
    return(z)},exo.sets,mc.cores = length(binsizes)))
  names(exo.densities) = binsizes

  pet.densities = suppressWarnings(mclapply(bins,function(b,seq.sets){
    z = list()
    z[[1]] = density.reads.per.strand.ratio(b,seq.sets[[1]])
    z[[2]] = density.reads.per.strand.ratio(b,seq.sets[[2]])
    return(z)},pet.sets,mc.cores = length(binsizes)))
  names(pet.densities) = binsizes
  
  exo.df = lapply(exo.densities,function(x){
    df1 = density2df(x[[1]])
    df2 = density2df(x[[2]])
    df1$Rep = 1
    df2$Rep = 2
    return(rbind(df1,df2))})
  names(exo.df) = binsizes

  pet.df = lapply(pet.densities,function(x){
    df1 = density2df(x[[1]])
    df2 = density2df(x[[2]])
    df1$Rep = 1
    df2$Rep = 2
    return(rbind(df1,df2))})
  names(pet.df) = binsizes

  exo.df = lapply(binsizes,function(b,seq.df){
    ss = seq.df[[as.character(b)]]
    ss$binSize = b
  return(ss)},exo.df)

  pet.df = lapply(binsizes,function(b,seq.df){
    ss = seq.df[[as.character(b)]]
    ss$binSize = b
  return(ss)},pet.df)

  exo.df = do.call(rbind,exo.df)
  exo.df$seq = "Exo"

  pet.df = do.call(rbind,pet.df) 
  pet.df$seq = "PET"
  
  df = rbind(exo.df,pet.df)
  df$seq = factor(df$seq)
  return(df)
}
  
density2df <- function(dens)return(data.frame("Fwd.Strand.Ratio" = dens$x ,density=dens$y)) 
  
plot.density <- function(binSize,exo.sets,pet.sets,genomeLength = seqlengths(exo.sets[[1]]))
{
  bins = create.bins(binSize,genomeLength)
  exo.densities = suppressWarnings(mclapply(exo.sets,function(x,bins)density.reads.per.strand.ratio(bins,x),bins,mc.cores = 2))
  pet.densities = suppressWarnings(mclapply(pet.sets,function(x,bins)density.reads.per.strand.ratio(bins,x),bins,mc.cores =2)) 
  
  exo.densities = lapply(1:length(exo.densities),function(j,dens){
    den = density2df(dens[[j]])
    den$Rep = j
    return(den)},exo.densities)
  exo.densities = do.call(rbind,exo.densities)
  exo.densities$seq = "ChIP-Exo"

  pet.densities = lapply(1:length(pet.densities),function(j,dens){
    den = density2df(dens[[j]])
    den$Rep = j
    return(den)},pet.densities)
  pet.densities = do.call(rbind,pet.densities)
  pet.densities$seq = "ChIP-Seq-PET"

  df = rbind(exo.densities,pet.densities)
  df$Rep = factor(df$Rep)
  df$seq = factor(df$seq)

  p <- ggplot(df,aes(Fwd.Strand.Ratio,density,colour = seq))+
    geom_line()+facet_grid(.~ Rep)+
    theme(legend.position = "top")
  return(p)
  
}
