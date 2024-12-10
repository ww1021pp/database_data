#BiocManager::install("recount3")
rm(list = ls())
setwd("/workspace/rsrch1/tmp/DataBase_datatable0618/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict")

library(recount3)
library(purrr)
library(dplyr)
library(ggplot2)

####download TCGA raw counts not run######
human_projects <- available_projects()

tcga_info = subset(
  human_projects,
  file_source == "tcga" & project_type == "data_sources"
)

head(tcga_info)


proj_info <- map(seq(nrow(tcga_info)), ~tcga_info[.x, ])

## create the RangedSummarizedExperiment. the create_rse function works on 
## one row a time 
rse_tcga <- map(proj_info, ~create_rse(.x))

saveRDS(rse_tcga,"TCGA_download.rawCounts.rds")

####################
matchedSample <- fread("../Sample_stastics/TumorMatchNATSampleInfo_withNATsampleGT25.txt")
matchedSample <- fread("../Sample_stastics/TumorMatchNATSampleInfo_withNATsampleGT25.txt")

#TCGA <- readRDS("TCGA_download.rawCounts.rds")

metaInfo=do.call(rbind,lapply(1:length(TCGA),function(x){
  data = as.data.frame(cbind(TCGA[[x]]@colData@rownames,do.call(cbind,TCGA[[x]]@colData@listData))) 
  
  return(data)  
})
)

matchedMeta <- metaInfo %>%
  select(V1,tcga.gdc_cases.project.project_id,tcga.cgc_sample_id,tcga.cgc_case_id,tcga.cgc_sample_sample_type,
         tcga.cgc_case_age_at_diagnosis,tcga.cgc_case_gender) %>% filter(tcga.cgc_sample_id %in% matchedSample$Sample.ID) %>% 
  group_by(tcga.cgc_sample_id) %>% dplyr::slice(1) %>% ungroup() %>% as.data.frame()

write.table(matchedMeta,"MatchSample.metaInfo.txt")
setdiff(matchedSample$Sample.ID,matchedMeta$tcga.cgc_sample_id)
# [1] "TCGA-A7-A0DC-01A" "TCGA-A7-A0DC-11A"


######rawcounts 
rawcounts_tmp <- lapply(1:length(TCGA),function(x){
  data <-as.data.frame(TCGA[[x]]@assays@data$raw_counts) %>% select(intersect(TCGA[[x]]@colData@rownames,matchedMeta$V1))
  return(data)
}
)
rawcounts <- do.call(cbind,rawcounts_tmp)

rownames(matchedMeta) = matchedMeta$V1
name_tmp = matchedMeta[colnames(rawcounts),"tcga.cgc_sample_id"]
names(rawcounts)=name_tmp
write.table(rawcounts,"matchedSample_rawCounts.txt",sep="\t",quote = F,row.names=F)

write.table(GeneName,"Gene_name.txt",row.names = F,quote = F)










