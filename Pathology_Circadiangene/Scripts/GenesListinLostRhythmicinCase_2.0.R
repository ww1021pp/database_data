rm(list=ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/Condition_CirGenesCorrelation/")

library(openxlsx)
library(data.table)
library(dplyr)
library(qdapTools)

df_table <- read.xlsx("./Data_res/LostRhyInCase/LostRhy_inCasecomparedControl.Genelist_table.xlsx")

comp_table <- read.xlsx("./Data_res/LostRhyInCase/LostRhymicInCase.GeneList_table.xlsx")

allCondList = as.list(df_table)

Bmalindex = grepl("B.*l.*KO",names(allCondList),ignore.case = T)
Cryindex = grepl("Cry.*KO",names(allCondList),ignore.case = T)
REVERB = grepl("REV",names(allCondList),ignore.case = T)
index = !(Bmalindex | Cryindex |REVERB)

NonKOList = allCondList[index]

ALF_binary = t(mtabulate(NonKOList)) %>% data.frame() %>% rownames_to_column(.,var = "Gene")

GeneinCondition=rowSums(ALF_binary[,-1])
ALF_binary$sum=GeneinCondition
test.list = lapply(unique(sort(ALF_binary$sum)),function(x){
  return(ALF_binary$Gene[which(ALF_binary$sum >= x)])
}
)

names(test.list)=paste0("CommonIngt",unique(sort(ALF_binary$sum)),"Condition")
data_list = list(
  LostRhymic_CaseConditionvsWT = ALF_binary,
  genelistinConds = list2DF(lapply(test.list, `length<-`, max(lengths(test.list))))
)
write.xlsx(data_list,"./Data_res/LostRhyInCase/NonCoreKO.LostRhymic_CaseConditionvsWT_2Sheets.xlsx")


Freq_con = as.data.frame(sapply(test.list,length)) %>% rownames_to_column(.,var="GeneinCondition")
names(Freq_con)[2]="Freq"
pdf("./plots_1219/LostRhythmic_plots/the histgram of frequency in NonKOLostCondition.pdf",height=4,width=4)
  ggplot(Freq_con, aes(x=GeneinCondition, y=Freq)) + 
  geom_bar(stat = "identity", width=0.3,fill="blue") +theme_classic()+
    xlab("Number of Conditions") +ylab("LostRhythmic Genes number in case") +coord_flip()
dev.off()
