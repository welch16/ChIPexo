

## this function only returns the characteristics of Landick's
## experiments for the different sequencing alternatives

edsn_tab <- function(what){
  stopifnot(what %in% c("exo","pet","set"))
  if(what == "exo"){
    edsn <- as.character(931:938)
    ip <- rep(c("Sig70","SigmaS"),4)
    condition <- rep(c("exp","stat"),each = 4)
    repl <- rep( rep(1:2,each =2),2)
    dt1 <- data.table(edsn,ip,condition,repl)
    edsn <- as.character(1310:1321)
    ip <- rep(c("Beta","Sig70","BetaPF"),4)
    condition <- rep( rep(c("rif0min","rif20min"),each = 3),2)
    repl <- rep(rep(1:2,each = 6))
    dt2 <- data.table(edsn,ip,condition,repl)
    dt <- rbind(dt1,dt2)    
  }else{
    edsn <- as.character(1396:1403)
    ip <- rep(c("Sig70","BetaPF"),4)
    condition <- rep(c("rif0min","rif20min"),each = 2, 2)
    repl <- rep(1:2,each = 4)
    dt <- data.table(edsn,ip,condition,repl)      
  }
  return(dt)
}
