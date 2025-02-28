rm(list=ls())
setwd("/workspace/rsrch1/tmp/DataBase_datatable0618/TCGA_NATSample_timePoints_predict")

library(data.table)
library(dplyr)
library(tibble)
sampleInfo <- read.delim("/workspace/rsrch1/tmp/TCGA/gdc_sample_sheet.2024-09-02.tsv")

#TCGA_count <- readRDS("/workspace/rsrch1/tmp/TCGA/count.rds")


cancer_summary <- data.frame(
  Cancer_Type = character(),
  Total_Samples = numeric(),
  Tumor_Samples = numeric(),
  Normal_Samples = numeric(),
  stringsAsFactors = FALSE
)



Normal_Samples <- sampleInfo %>% 
                  filter(Sample.Type == "Solid Tissue Normal") %>%
                  group_by(Project.ID) %>% summarise(Normal_sample=n())

Tumor_Samples <- sampleInfo %>% 
  filter(Sample.Type != "Solid Tissue Normal") %>%
  group_by(Project.ID) %>% summarise(Tumor_sample=n()) %>% left_join(.,y=Normal_Samples)


Normal_Samples_info <- sampleInfo %>% filter(Sample.Type == "Solid Tissue Normal")

matchedCase_sample_info <- sampleInfo %>% filter(Case.ID %in% Normal_Samples_info$Case.ID & 
                                                 Sample.Type != "Solid Tissue Normal") %>%
                           group_by(Case.ID) %>% arrange(Sample.ID) %>% slice(1) %>% ungroup()
matchedNAT_sample_info <-  Normal_Samples_info %>% filter(Case.ID %in% intersect(Normal_Samples_info$Case.ID,matchedCase_sample_info$Case.ID))

TumorNAT_matchedSampleInfo <- rbind(matchedCase_sample_info,matchedNAT_sample_info)

tmp <- TumorNAT_matchedSampleInfo %>% group_by(Case.ID) %>% summarise(n())

matchedNormal_Samples <- TumorNAT_matchedSampleInfo %>% 
  filter(Sample.Type == "Solid Tissue Normal") %>%
  group_by(Project.ID) %>% summarise(Normal_sample=n())

matchedTumor_Samples <- TumorNAT_matchedSampleInfo %>% 
  filter(Sample.Type != "Solid Tissue Normal") %>%
  group_by(Project.ID) %>% summarise(Tumor_sample=n()) %>% left_join(.,y=matchedNormal_Samples)

head(matchedTumor_Samples);
head(Tumor_Samples)

#Tumor_Samples[is.na(Tumor_Samples)]=0
#matchedTumor_Samples[is.na(matchedTumor_Samples)]=0

longallSamples <- reshape2::melt(Tumor_Samples, id.vars = c("Project.ID"), variable.name = "SampleType")
longMatchedSamples <- reshape2::melt(matchedTumor_Samples,id.vars = c("Project.ID"), variable.name = "SampleType")

library(ggplot2)

# Stacked
pdf("./Sample_stastics/all_sampleStastics.pdf",height = 6,width = 4)
ggplot(longallSamples, aes(fill=SampleType, y=value, x=Project.ID)) + 
  geom_bar(position="stack", stat="identity") + theme_classic()+
  geom_text(aes(label = value), size = 2, position = "stack")+
  coord_flip() + theme(legend.position = "top",
                       axis.text.y=element_text(hjust=0),
                       axis.text = element_text(colour = "black",size=8))
dev.off()



pdf("./Sample_stastics/matchedTumor_NATsampleStastics.pdf",height = 6,width = 4)
ggplot(longMatchedSamples, aes(fill=SampleType, y=value, x=Project.ID)) + 
  geom_bar(position="stack", stat="identity") + theme_classic()+
  geom_text(aes(label = value), size = 2, position = "stack")+
  coord_flip() + theme(legend.position = "top",
                       axis.text.y=element_text(hjust=0),
                       axis.text = element_text(colour = "black",size=8))
dev.off()


table(matchedTumor_Samples$Normal_sample >=30)

ProjID <- matchedTumor_Samples[which(matchedTumor_Samples$Normal_sample>=25),]

write.table(TumorNAT_matchedSampleInfo,
            "./Sample_stastics/TumorMatchNATSampleInfo.txt",sep="\t",quote = F,row.names = F)

write.table(TumorNAT_matchedSampleInfo[which(TumorNAT_matchedSampleInfo$Project.ID %in% ProjID$Project.ID),],
            "./Sample_stastics/TumorMatchNATSampleInfo_withNATsampleGT25.txt",sep="\t",quote = F,row.names = F)









