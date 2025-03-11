rm(list=ls())
setwd("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/")

library(openxlsx)
library(dplyr)
library(data.table)

genelist = read.xlsx("./LostRhyInCase/NonCoreKO.LostRhymic_CaseConditionvsWT_2Sheets.xlsx",sheet = 2)

bedfile = read.delim("/workspace/rsrch2/common_data/Refgenome/mm10/GenCode.vM25.mm10.Gene_bed/GenCodeVM25.gene.tssup2.5kbDn1kb.padded.filtered.bed",header = F)
names(bedfile)=c("Chr","Start","End","Gene","Score","Strand")


bed_res <- lapply(8:dim(genelist)[2],function(x){
  genes=na.omit(genelist[,x])
  befile_file = bedfile %>% filter(Gene %in% genes)
  filename = paste0("./LostRhyInCase/Lost_TF.BedFiles/",names(genelist)[x],".bed")
  #filename=paste0("./GainRhyInCase/Gain_TF.BedFiles/CommonInGT",names(genelist)[x],"condition.bed")
  write.table(befile_file,file = filename,sep = "\t",row.names = F,col.names = F,quote = F)
  })
