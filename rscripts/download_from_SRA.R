
rm(list = ls())

library(SRAdb)
library(DBI)

srafile <- getSRAdbFile()
con <- dbConnect(RSQLite::SQLite(), srafile)

listSRAfile('SRR770743', con)
listSRAfile('SRR770744', con)
listSRAfile('SRA046523',con)


listSRAfile('SRP038893',con)

getSRAfile('SRR770743',con,fileType='sra')
getSRAfile('SRR770744',con,fileType='sra')
getSRAfile('SRA046523',con,fileType='sra')


