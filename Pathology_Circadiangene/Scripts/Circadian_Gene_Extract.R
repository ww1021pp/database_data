m(list=ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/Pathology_ChangedGenes/")

library(openxlsx)
library(data.table)
library(dplyr)
library(psych)
library(qdapTools)
################RNA-seq data###########
cond_df.Liver_raw <- read.xlsx("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/MetacycleR/RNA-seq.condition.xlsx",sheet = 2)
all_con <- read.xlsx("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/MetacycleR/RNA-seq.condition.xlsx",sheet = 1)
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
  if(all_con[which(all_con$gseID==gsub(".*(GSE[0-9]+).*","\\1",i))[1],"LibraryType"]=="Microarray"){pvalue_col="JTK_pvalue"}
  
  amp_col=colnames(tmp)[grep(gsub("pvalue","amp",pvalue_col),ignore.case = T,colnames(tmp))]
  phase=colnames(tmp)[grep(paste0(gsub("_pvalue","",pvalue_col),".*phase"),ignore.case = T,colnames(tmp))]
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

Pathology.meta = c("GSE132103_metaout_CHRONIC_EtOH.txt",
                 "GSE133989_metaout_H2O_RHYTHM.txt",
                 "GSE52333_MicroArrayHFD/GSE52333_metaout.Liver.HFD.txt",
                 "GSE73222_LungCancer/GSE73222_metaout.Liver.TB.txt",
                 "GSE118967_metaout_AR.txt",
                 "GSE93903_MicroArrayAge/GSE93903_metaout.liver.OldND.txt",
                 "GSE143528_metaout_RPF_WT.txt",
                 "GSE220864_metaout_KO_LI.txt",
                 "PRJEB16726_metaout_Abx.txt"
                 ) 
df_Pathology = all_con %>% filter(metaoutfile %in% Pathology.meta) %>% mutate(
  pathfile = case_when(LibraryType =="RNA-Seq" ~ paste0(metafile_path,"/",gseID,"/",metaoutfile),
                       LibraryType =="Microarray" ~ paste0(Arraypath,"/",metaoutfile)
  )
  )


TRF.meta =c("GSE158600_TRFWTKO/GSE158600_metaout.TRF_WT.txt",
            "GSE135898_metaout_Bmal1_WT_RF.txt",
            "GSE73552_metaout_RNA_WT_RF.txt",
            "GSE93903_MicroArrayAge/GSE93903_metaout.liver.YoungCR.txt",
            "GSE93903_MicroArrayAge/GSE93903_metaout.liver.OldCR.txt")

df_TRF = all_con %>% filter(metaoutfile %in% TRF.meta) %>% mutate(
  pathfile = case_when(LibraryType =="RNA-Seq" ~paste0(metafile_path,"/",gseID,"/",metaoutfile),
                       LibraryType =="Microarray" ~paste0(Arraypath,"/",metaoutfile)
  )
)

Pathology=list()
list_tmp=lapply(df_Pathology$pathfile[2:8],circadian_list)
GSE220864 = fread(df_Pathology$pathfile[9]) %>% filter(JTK_pvalue <= 0.05) %>%
  mutate(Gene=Gene,pvalue=JTK_pvalue,amp=JTK_amplitude) %>% group_by(Gene) %>% arrange(pvalue) %>%
   dplyr::slice(1)%>% ungroup() %>% dplyr::select(Gene,pvalue,amp)
Pathology[2:8]=list_tmp
Pathology[[9]]=GSE220864
names(Pathology)[2:9]=df_Pathology$metaoutfile[2:9]

PRJEB16726 = fread(df_Pathology$pathfile[1]) %>% filter(JTK_pvalue <= 0.05) %>%
  mutate(Gene=Gene,pvalue=JTK_pvalue,amp=JTK_amplitude) %>% group_by(Gene) %>% arrange(pvalue) %>% 
  dplyr::slice(1)%>% ungroup()%>% dplyr::select(Gene,pvalue,amp)
Pathology[[1]]=PRJEB16726
names(Pathology)[1]=df_Pathology$metaoutfile[1]

Pathology_con = lapply(1:length(Pathology),function(x){
  tmp=Pathology[[x]]
  con = gsub(".*metaout_(.*)\\.txt","\\1",names(Pathology)[x])
  gseID=gsub(".*(GSE[0-9]+).*","\\1",names(Pathology)[x])
 tmp$gseID=gseID
 tmp$con=con
 tmp$Group="Pathology"
 return(tmp)
})
names(Pathology_con)=df_Pathology$metaoutfile

saveRDS(Pathology_con,"NineCondition_CircadianList.rds")


################TRF circadian genes####

TRF_cir=lapply(df_TRF$pathfile,circadian_list)
names(TRF_cir)=df_TRF$metaoutfile
TRF_con = lapply(1:length(TRF_cir),function(x){
  tmp=TRF_cir[[x]]
  con = gsub(".*metaout_(.*)\\.txt","\\1",names(TRF_cir)[x])
  gseID=gsub(".*(GSE[0-9]+).*","\\1",names(TRF_cir)[x])
  tmp$gseID=gseID
  tmp$con=con
  tmp$Group="NRF"
  return(tmp)
})
names(TRF_con)=df_TRF$metaoutfile
saveRDS(TRF_con,"FiveNRFCircadianList.rds")



