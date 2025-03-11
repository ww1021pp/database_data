rm(list=ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/")

library(openxlsx)
library(data.table)
library(dplyr)

metaInfo = read.delim("Cistrom_mm10_tranfac_QC.txt",header = T) %>% select(2,3,5,6)

path="~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/LostRhyInCase/Lost_TF.BedFiles/GIGGLE_res/"
#path="~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/Gain_TF.BedFiles/Gainrhy_GIGGLE_res/"
output=path
data <- list.files(path=path,pattern = "csv",full.names = T)
col_name=c("file","file_size","overlaps","odds_ratio","fishers_two_tail","fishers_left_tail","fishers_right_tail","GIGGLE_score")

df_list <- lapply(data,function(x){
  print(x)
  df = as.data.frame(read.csv(x,sep = "\t",header = F) %>% select(1:8))
  df_raw=df[-1,] 
  names(df_raw)=col_name
  df_res = df_raw %>% mutate(External_id=gsub(".*\\/(.*)_sort.*.gz","\\1",file),
                             GIGGLE_score = as.numeric(GIGGLE_score)) %>% left_join(.,y=metaInfo)
  return(df_res)
  })

names(df_list)= gsub("_GIGGLE_res_Top1K.csv","",basename(data))

df_res <- lapply(df_list,function(x){
  tmp = x
  res_tmp = tmp %>% filter(grepl("[Ll]iver",Ontology)) %>% group_by(Factor) %>% 
    arrange(desc(as.numeric(GIGGLE_score))) %>% dplyr::slice(1) %>% ungroup() %>% arrange(desc(GIGGLE_score))
  return(res_tmp)
  
})

names(df_res)=names(df_list)
write.xlsx(df_res,paste0(output,"GIGGLE_res.xlsx"))


