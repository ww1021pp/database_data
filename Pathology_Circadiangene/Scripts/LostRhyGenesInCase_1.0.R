rm(list=ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/Condition_CirGenesCorrelation/")

library(openxlsx)
library(data.table)
library(dplyr)
library(psych)
library(qdapTools)
################RNA-seq data###########
cond_df.Liver_raw <- read.xlsx("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/MetacycleR/RNA-seq.condition.xlsx",sheet = 2)
temp <- unique(cond_df.Liver_raw$GseID[duplicated(cond_df.Liver_raw$GseID)])
index=which(temp == "GSE155160")
temp=temp[-index]
metafile_path = "/workspace/rsrch1/tmp/DataBase_datatable0618/MetaCycle/APA_results/MetacycleR/metaCycle_Res/MetaResults_Expression"
meta.results = list.files(path=metafile_path,
                          pattern = "meta2d",recursive = T)
df_tmp1 = data.frame(
  GseID=gsub("(.*)\\/.*\\/.*","\\1",meta.results),
  folder = gsub(".*\\/(.*)\\/.*","\\1",meta.results),
  ExpMeta = paste0(metafile_path,"/",meta.results),
  library = "RNA-seq"
)

############Microarrary data for liver#### ~/project/CirDatabase_020824/Kristin.SeriesePaper
Arraypath="~/project/CirDatabase_020824/Kristin.SeriesePaper"
metaout = list.files(path = Arraypath,pattern = "metaout",recursive = T)
df_tmp2 = data.frame(
  GseID = gsub("(.*)_.*\\/.*","\\1",metaout),
  folder = gsub(".*\\/(.*)","\\1",metaout),
  ExpMeta = paste0(Arraypath,"/",metaout),
  library = "Microarray"
)

df_tmp=rbind(df_tmp1,df_tmp2)
cond_df.Liver_raw$FolderName = gsub(" ","",cond_df.Liver_raw$FolderName)
cond_df.Liver = inner_join(cond_df.Liver_raw, df_tmp,by=c("GseID"="GseID","FolderName"="folder")) %>% filter(GseID %in% temp)

Ensem_Gene <- read.delim("/workspace/rsrch1/tmp/DataBase_datatable0618/metaInfo_Folder/3Organism_geneID.txt") 
EnsemGene=Ensem_Gene %>% dplyr::filter(Organism == "Mus musculus") %>% dplyr::select(ensembl_transcript_id,Gene) %>% dplyr::filter(!is.na(Gene)) %>% unique()
rownames(EnsemGene)=EnsemGene$ensembl_transcript_id

circadian_list <- function(i){
  tmp <- fread(i)
  pvalue_col=intersect(c("meta2d_pvalue","JTK_pvalue","LS_pvalue"),colnames(tmp))[1]
  colnames(tmp)[1]="EnsembleID"
  pvalue=ifelse(pvalue_col=="LS_pvalue",0.1,0.05)
  amp_col=colnames(tmp)[grep(gsub("pvalue","amp",pvalue_col[1]),ignore.case = T,colnames(tmp))]
  phase=colnames(tmp)[grep(paste0(gsub("_pvalue","",pvalue_col[1]),".*phase"),ignore.case = T,colnames(tmp))]
  if(grepl("^ENSMUST",tmp$EnsembleID[1])){
    res <- tmp %>% filter( .[[pvalue_col]] < pvalue)  %>% mutate(Gene=EnsemGene[EnsembleID,"Gene"],
                                                                                              pvalue=.[[pvalue_col]],
                                                                                              amp=.[[amp_col]]
                                                                                              ) %>% filter(Gene !="") %>% arrange(Gene,pvalue) %>% group_by(Gene) %>% dplyr::slice(1) %>% ungroup %>% dplyr::select(Gene,pvalue,amp) %>% as.data.frame()
    } else{
    res <- tmp %>% filter( .[[pvalue_col]] < pvalue) %>% mutate(Gene=gsub("\\s.*","",EnsembleID),
                                                                pvalue=.[[pvalue_col]],
                                                                amp=.[[amp_col]]) %>% filter(Gene !="")  %>% 
      filter(Gene !="") %>% arrange(Gene,pvalue) %>% group_by(Gene) %>% dplyr::slice(1) %>% ungroup %>% dplyr::select(Gene,pvalue,amp) %>% as.data.frame()
    }
  
  rownames(res)=res$Gene
  return(res)
}

compareWT = function(x){
  df <- cond_df.Liver %>% filter(GseID == x) %>% arrange(Subgroup)
  res <- lapply(df$ExpMeta,circadian_list)
  names(res)=df$Subgroup
  
  group_num=nrow(df)/2
  
  if(df$Subgroup[nrow(df)] == "Ctrl"){
    test=do.call('rbind',lapply(1:(nrow(df)-1),function(j){
      gainRhy = data.frame(
        GseID=x,     
        dataset=paste0(x,"_",df$Condition[j],".",df$Diet[j],"VS",df$Condition[nrow(df)],".",df$Diet[nrow(df)]),
        gain = setdiff(res[[nrow(df)]]$Gene,res[[j]]$Gene),
        pvalue=res[[nrow(df)]][setdiff(res[[nrow(df)]]$Gene,res[[j]]$Gene),"pvalue"],
        amp = res[[nrow(df)]][setdiff(res[[nrow(df)]]$Gene,res[[j]]$Gene),"amp"]
    )
        
      return(gainRhy)
    }
    )
    )
    return(test)
  }
  if(grepl("Ctrl1",df$Subgroup[group_num+1])){
    test=do.call('rbind',lapply(1:group_num,function(j){
      gainRhy = data.frame(
        GseID=x,     
        dataset=paste0(x,"_",df$Condition[j],".",df$Diet[j],"VS",df$Condition[group_num+j],".",df$Diet[group_num+j]),
        gain = setdiff(res[[group_num+j]]$Gene,res[[j]]$Gene),
        pvalue=res[[group_num+j]][setdiff(res[[group_num+j]]$Gene,res[[j]]$Gene),"pvalue"],
        amp = res[[group_num+j]][setdiff(res[[group_num+j]]$Gene,res[[j]]$Gene),"amp"]
      )
      return(gainRhy)
    }
    )
    )
    return(test)
  }
  else{next}
}




Res_setDiff1_9 <- lapply(temp[1:9],compareWT)
names(Res_setDiff1_9)=temp[1:9]

Res_setDiff11_19<- lapply(temp[c(11:19)],compareWT)
names(Res_setDiff11_19)=temp[c(11:19)]


##########GSE143528 TRFvsWT####
data <- cond_df.Liver %>% filter(GseID == "GSE143528") %>% arrange(Subgroup)
res_data <- lapply(data$ExpMeta,circadian_list)
names(res_data)=data$Subgroup
gain_temp = data.frame(
  GseID="GSE143528",     
  dataset=paste0("GSE143528","_","TRFvsWT"),
  gain = setdiff(res_data[[3]]$Gene,res_data[[4]]$Gene),
  pvalue=res_data[[3]][setdiff(res_data[[3]]$Gene,res_data[[4]]$Gene),"pvalue"],
  amp = res_data[[3]][setdiff(res_data[[3]]$Gene,res_data[[4]]$Gene),"amp"]
  
)
Res_setDiff11_19$GSE143528=rbind(Res_setDiff11_19$GSE143528,gain_temp)



###################GSE158600###########
data <- cond_df.Liver %>% filter(GseID == temp[20]) %>% arrange(Subgroup)
res_data <- lapply(data$ExpMeta,circadian_list)
names(res_data)=data$Subgroup
test=do.call('rbind',lapply(1:2,function(j){
  gainRhy = data.frame(
    GseID=temp[20],     
    dataset=paste0(temp[20],"_",data$Condition[j],".",data$Diet[j],"VS",data$Condition[3],".",data$Diet[3]),
    gain = setdiff(res_data[[3]]$Gene,res_data[[j]]$Gene),
    pvalue=res_data[[3]][setdiff(res_data[[3]]$Gene,res_data[[j]]$Gene),"pvalue"],
    amp = res_data[[3]][setdiff(res_data[[3]]$Gene,res_data[[j]]$Gene),"amp"]
    )
  return(gainRhy)
}
)
)

test_NRF=do.call('rbind',lapply(5:6,function(j){
  gainRhy = data.frame(
    GseID=temp[20],     
    dataset=paste0(temp[20],"_",data$Condition[j],".",data$Diet[j],"VS",data$Condition[4],".",data$Diet[4]),
    gain = setdiff(res_data[[4]]$Gene,res_data[[j]]$Gene),
    pvalue=res_data[[4]][setdiff(res_data[[4]]$Gene,res_data[[j]]$Gene),"pvalue"],
    amp = res_data[[4]][setdiff(res_data[[4]]$Gene,res_data[[j]]$Gene),"amp"]
  )
  return(gainRhy)
}
)
)

BmalTRFvsAL = data.frame(
  GseID=temp[20],     
  dataset=paste0(temp[20],"_",data$Condition[4],".",data$Diet[4],"VS",data$Condition[3],".",data$Diet[3]),
  gain = setdiff(res_data$Ctrl$Gene,res_data$NRF_Ctrl$Gene),
  pvalue=res_data$Ctrl[setdiff(res_data$Ctrl$Gene,res_data$NRF_Ctrl$Gene),"pvalue"],
  amp = res_data$Ctrl[setdiff(res_data$Ctrl$Gene,res_data$NRF_Ctrl$Gene),"amp"]
)
Res_setDiff20 = list(rbind(test,test_NRF,BmalTRFvsAL))
names(Res_setDiff20)=temp[20]


Res_setDiff21_26=lapply(temp[21:26],compareWT)
names(Res_setDiff21_26)=temp[21:26]


#######temp[10] GSE135898#######

data <- cond_df.Liver %>% filter(GseID == "GSE135898") %>% arrange(Subgroup)
res <- lapply(data$ExpMeta,circadian_list)
names(res)=data$Subgroup

Bmal=do.call('rbind',lapply(1:2,function(j){
  gainRhy = data.frame(
    GseID=temp[10],     
    dataset=paste0(temp[10],"_",data$Condition[j],".",data$Diet[j],"VS",data$Condition[2+j],".",data$Diet[2+j]),
    gain = setdiff(res[[2+j]]$Gene,res[[j]]$Gene),
    pvalue=res[[2+j]][setdiff(res[[2+j]]$Gene,res[[j]]$Gene),"pvalue"],
    amp = res[[2+j]][setdiff(res[[2+j]]$Gene,res[[j]]$Gene),"amp"]
    
  )
  return(gainRhy)
}
)
)

Cry=do.call('rbind',lapply(5:6,function(j){
  gainRhy = data.frame(
    GseID=temp[10],     
    dataset=paste0(temp[10],"_",data$Condition[j],".",data$Diet[j],"VS",data$Condition[2+j],".",data$Diet[2+j]),
    gain = setdiff(res[[2+j]]$Gene,res[[j]]$Gene),
    pvalue=res[[2+j]][setdiff(res[[2+j]]$Gene,res[[j]]$Gene),"pvalue"],
    amp = res[[2+j]][setdiff(res[[2+j]]$Gene,res[[j]]$Gene),"amp"]
  )
  return(gainRhy)
}
)
)

Bmal_TRFvsAL = data.frame(
  GseID=temp[10],
  dataset=c(rep(paste0(temp[10],"_BmalTRFVSAL"),length(setdiff(res$BmalCtrl1$Gene,res$BmalCtrl2$Gene))),
            rep(paste0(temp[10],"_CryTRFVSAL"),length(setdiff(res$CryCtrl1$Gene,res$CryCtrl2$Gene)))),
  gain=c(setdiff(res$BmalCtrl1$Gene,res$BmalCtrl2$Gene),setdiff(res$CryCtrl1$Gene,res$CryCtrl2$Gene)),
  pvalue=c(res$BmalCtrl1[setdiff(res$BmalCtrl1$Gene,res$BmalCtrl2$Gene),"pvalue"],res$CryCtrl1[setdiff(res$CryCtrl1$Gene,res$CryCtrl2$Gene),"pvalue"]),
  amp = c(res$BmalCtrl1[setdiff(res$BmalCtrl1$Gene,res$BmalCtrl2$Gene),"amp"],res$CryCtrl1[setdiff(res$CryCtrl1$Gene,res$CryCtrl2$Gene),"amp"])
)


Res_setDiff10 = list(rbind(Bmal,Cry,Bmal_TRFvsAL))
names(Res_setDiff10)=temp[10]

#####################Old vs Young GSE93903###
data <- cond_df.Liver %>% filter(GseID == "GSE93903") %>% arrange(Subgroup)
res <- lapply(data$ExpMeta,circadian_list)
names(res)=data$Subgroup
OldvsYoung = data.frame(
  GseID="GSE93903",     
  dataset=paste0("GSE93903","_",data$Condition[4],".",data$Diet[4],"VS",data$Condition[3],".",data$Diet[3]),
  gain = setdiff(res[[3]]$Gene,res[[4]]$Gene),
  pvalue=res[[3]][setdiff(res[[3]]$Gene,res[[4]]$Gene),"pvalue"],
  amp = res[[3]][setdiff(res[[3]]$Gene,res[[4]]$Gene),"amp"]
)

Res_setDiff21_26$GSE93903=rbind(Res_setDiff21_26$GSE93903,OldvsYoung)

####################add new datasets GSE220864#########3
data <- cond_df.Liver %>% filter(GseID == "GSE220864") %>% arrange(Subgroup)
Ctrl <- fread("/workspace/rsrch1/tmp/DataBase_datatable0618/MetaCycle/APA_results/MetacycleR/metaCycle_Res/MetaResults_Expression/GSE220864/GSE220864_metaout_WT_LI.txt",header = T)
Case <- fread("/workspace/rsrch1/tmp/DataBase_datatable0618/MetaCycle/APA_results/MetacycleR/metaCycle_Res/MetaResults_Expression/GSE220864/GSE220864_metaout_KO_LI.txt",header = T)
CirList=list(Ctrl = Ctrl %>% filter(JTK_pvalue <= 0.05 & Gene !="-") %>% group_by(Gene) %>% arrange(JTK_pvalue) %>% dplyr::slice(1) %>% ungroup() %>% dplyr::select(Gene,JTK_pvalue,JTK_amplitude) %>% as.data.frame(),
             Case = Case%>% filter(JTK_pvalue <= 0.05 & Gene !="-") %>% group_by(Gene) %>% arrange(JTK_pvalue) %>% dplyr::slice(1) %>% ungroup() %>%  dplyr::select(Gene,JTK_pvalue,JTK_amplitude) %>% as.data.frame())
rownames(CirList$Ctrl)=CirList$Ctrl$Gene
rownames(CirList$Case)=CirList$Case$Gene
LostRhy=data.frame(
  GseID="GSE220864",     
  dataset=paste0("GSE220864","_",data$Condition[1],".",data$Diet[1],"VS",data$Condition[2],".",data$Diet[2]),
  LostRhy= setdiff(CirList$Ctrl$Gene,CirList$Case$Gene),
  pvalue=CirList$Ctrl[setdiff(CirList$Ctrl$Gene,CirList$Case$Gene),"JTK_pvalue"],
  amp = CirList$Ctrl[setdiff(CirList$Ctrl$Gene,CirList$Case$Gene),"JTK_amplitude"]
)
write.table(LostRhy,"GSE220864_CaseLostRhy_022025.txt",sep = "\t",quote = F,)  



GainRhy=data.frame(
  GseID="GSE220864",     
  dataset=paste0("GSE220864","_",data$Condition[1],".",data$Diet[1],"VS",data$Condition[2],".",data$Diet[2]),
  GainRhy= setdiff(CirList$Case$Gene,CirList$Ctrl$Gene),
  pvalue=CirList$Case[setdiff(CirList$Case$Gene,CirList$Ctrl$Gene),"JTK_pvalue"],
  amp = CirList$Case[setdiff(CirList$Case$Gene,CirList$Ctrl$Gene),"JTK_amplitude"]
)
write.table(GainRhy,"GSE220864_CaseGainRhy_022025.txt",sep = "\t",quote = F,)  



####################all liver condition vs Corresponding WT###

all_list = c(Res_setDiff1_9,Res_setDiff10,Res_setDiff11_19,Res_setDiff20,Res_setDiff21_26)
saveRDS(all_list,"./Data_res/LostRhyInCase/LostRhythmic_CasecomparedwithCtrlP0.05amppvalue.rds")

all_list=readRDS("./Data_res/LostRhyInCase/LostRhythmic_CasecomparedwithCtrlP0.05amppvalue.rds")
all_cond <- do.call(rbind,all_list)
names(all_cond)[3]="LostRhy"
write.table(all_cond,"./Data_res/LostRhyInCase/LostRhy_rhyGenesinCasecomparedWT_ampPvalue.txt",sep="\t",quote = F)
