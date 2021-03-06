

## this function only returns the characteristics of Landick's
## experiments for the different sequencing alternatives

edsn_tab <- function(what){
  stopifnot(what %in% c("exo","pet","set"))
  if(what == "exo"){
    edsn <- as.character(1310:1321)
    ip <- rep(c("Beta","Sig70","BetaPF"),4)
    condition <- rep( rep(c("rif0min","rif20min"),each = 3),2)
    repl <- rep(rep(1:2,each = 6))
    dt <- data.table(edsn,ip,condition,repl)
  }else{
    edsn <- as.character(1396:1403)
    ip <- rep(c("Sig70","BetaPF"),4)
    condition <- rep(c("rif0min","rif20min"),each = 2, 2)
    repl <- rep(1:2,each = 4)
    dt <- data.table(edsn,ip,condition,repl)      
  }
  return(dt)
}


edsn_tab_old <- function(what){
  stopifnot(what %in% c("exo","pet","set"))
  if(what == "exo"){
    ## exo
#    edsn <- as.character(c("931"),"933"))
    edsn <- "931"
    ip <- rep("Sig70",1)
    condition <- rep("aerobic",1)
    growth <- rep("exp",1)
    repl <- 1
    dt <- data.table(edsn,ip,condition,growth,repl)    
  }else if(what == "pet"){    
    ## pet
    edsn <- "788"
    ip <- "Sig70"
    condition <- "aerobic"
    growth <- "exp"
    repl <- 1
    dt <- data.table(edsn,ip,condition,growth,repl)    
  }else{
    ## set
    edsn <- "80"
    ip <- "Sig70"
    condition <- "aerobic"
    growth <- "exp"
    repl <- 1
    dt <- data.table(edsn,ip,condition,growth,repl)        
  }
  return(dt)

}

## function that returns the right read files given a directory

filter_reads <- function(dir,char,isPET)
{
  files <- list.files(dir)
  files <- files[grep(char[,(edsn)],files)]
  files <- files[grep("bai",files,invert = TRUE)]  
  if(isPET){
    files <- files[grep("filter",files)]
  }
  files <- files[grep("run",files,invert = TRUE)]
  return(file.path(dir,files))
}
