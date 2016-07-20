
rm(list = ls())

library(mosaics)
library(parallel)
library(dpeak)
library(data.table)
library(GenomicRanges)
library(GenomicAlignments)

indir = "/p/keles/ChIPexo/volume6/K12/subsample_resolution_sensitivity"

files = list.files(indir,full.names = TRUE,recursive = TRUE)
files = files[grep("bins",files,invert = TRUE)]
files = files[grep("bin",files,invert = TRUE)]
files = files[grep(".pl",files,invert = TRUE,fixed = TRUE)]

reads = files[grep("subsample_1",files)]

mc = 24

peaks = files[grep("RData",files)]
peaks = peaks[grep("peaks",peaks)]

load(peaks) ## peaks

rm(files)

dr = tempdir()
files = vapply(seq_len(length(unlist(peaks))),
  function(x)tempfile(pattern = "peak",tmpdir = dr,fileext = ".bed"),"")

sqnms = rep(c("gi|49175990|ref|NC_000913.2|","U00096"),each = 4)

generate_bed = function(peak,file,sqnms = "U00096")
{
  mosaics::export(peak,filename = file,type = "bed")
  tab = fread(file)
  tab[,V1 := sqnms]
  tab[,V2 := V2 + 1]
  tab[,V3 := V3 + 1]
  write.table(tab,file = file,quote = FALSE,
              row.names = FALSE,col.names = FALSE)
  tab
}
   
tabs = mapply(generate_bed,unlist(peaks),files,sqnms,SIMPLIFY = FALSE)

## ChIP-exo
exoidx = 1:4
petidx = 5:6
setidx = 7:8
method = "separate"
nTop = 100

exo = mapply(dpeakRead,peakfile = files[exoidx],readfile = reads[exoidx],
  fileFormat = "sam",fragLen = 150, PET = FALSE,parallel = TRUE,
  nCore = mc,SIMPLIFY = FALSE)

exofit = lapply(exo,dpeakFit,estDeltaSigma = method,lbDelta = 5, lbSigma = 5,
  nCore = mc,maxComp = 5)

exofiles = vapply(seq_along(exofit),
  function(x)file.path(indir,"ChIPexo","sites",paste0("TFBS_",x,".txt")),"")

mapply(export,exofit,type = "txt",filename = exofiles,SIMPLIFY = FALSE)

## ChIP-seq PET
pet = mapply(dpeakRead,peakfile = files[petidx],readfile = reads[petidx],
  fileFormat = "eland_result",PET = TRUE,parallel = TRUE,
  nCore = mc,SIMPLIFY = FALSE)

petfit = lapply(pet,dpeakFit,estDeltaSigma = method,lbDelta = 5, lbSigma = 5,
  nCore = mc,maxComp = 5)

petfiles = vapply(seq_along(petfit),
  function(x)file.path(indir,"ChIPseq_PET","sites",paste0("TFBS_",x,".txt")),"")

mapply(export,petfit,type = "txt",petfiles,SIMPLIFY = FALSE)


## ChIP-seq SET
set = mapply(dpeakRead,peakfile = files[setidx],readfile = reads[setidx],
  fileFormat = "eland_result",PET = FALSE,parallel = TRUE,
  nCore = mc,fragLen = 150,
  SIMPLIFY = FALSE)

setfit = lapply(set,dpeakFit,estDeltaSigma = method,lbDelta = 5, lbSigma = 5,
  nCore = mc,maxComp = 5)

setfiles = vapply(seq_along(setfit),
  function(x)file.path(indir,"ChIPseq_SET","sites",paste0("TFBS_",x,".txt")),"")

mapply(export,setfit,type = "txt",setfiles,SIMPLIFY = FALSE)

save(exofit,file = file.path(indir,"ChIPexo","sites",paste0("TFBS_",method,
              ".RData")))

save(petfit,file = file.path(indir,"ChIPseq_PET","sites",
              paste0("TFBS_",method,
              ".RData")))

save(setfit,file = file.path(indir,"ChIPseq_SET","sites",
              paste0("TFBS_",method,
              ".RData")))
