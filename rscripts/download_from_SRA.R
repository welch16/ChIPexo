
rm(list = ls())

library(SRAdb)
library(DBI)

srafile <- getSRAdbFile()
con <- dbConnect(RSQLite::SQLite(), srafile)





listSRAfile("SRX248184",con)
listSRAfile("SRX248185",con)


listSRAfile('SRR770743', con)
listSRAfile('SRR770744', con)
listSRAfile('SRA046523',con)


listSRAfile('SRP038893',con)

getSRAfile('SRR770743',con,fileType='sra')
getSRAfile('SRR770744',con,fileType='sra')
getSRAfile('SRA046523',con,fileType='sra')

getSRAfile("SRX248184",con,fileType = "sra")
getSRAfile("SRX248185",con,fileType = "sra")


listSRAfile("SRA067908",con)$run

library(parallel)

mclapply(listSRAfile("SRA067908",con)$run,
       function(x)getSRAfile(x,con,fileType = "sra"),mc.cores =detectCores())
