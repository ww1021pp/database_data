rm(list=ls())
setwd("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/LostRhyInCase/LostRhy_commonTCGA")
library(openxlsx)
library(tibble)
library(dplyr)



###1305 1366 1314

Ensem_gene <- read.delim("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/Ensemble_gene_mappingHsp.txt") %>% 
  group_by(SYMBOL) %>% dplyr::slice(1) %>% ungroup() %>% as.data.frame()
rownames(Ensem_gene)=Ensem_gene$ENSEMBL

Mouse_LostrhyCase <- read.xlsx("../NonCoreKO.LostRhymic_CaseConditionvsWT_2Sheets.xlsx",sheet = 2)%>% as.list(na.omit(.))

TCGA = get(load("~/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/OUT/NATvsTumore.sample.FiT.OUT_Status.RData"))
Liver_NAT=TCGA$`TCGA-LIHC-NAT`$data.fit
Liver_Tumor = TCGA$`TCGA-LIHC-Tumor`$data.fit


Tumor_rhyGenes=Liver_Tumor %>% filter(pval < 0.05) 
NAT_rhyGenes = Liver_NAT %>% filter(pval < 0.05)

LostInTumor = na.omit(Ensem_gene[setdiff(rownames(NAT_rhyGenes),rownames(Tumor_rhyGenes)),"SYMBOL"])
write.table(LostInTumor,"LostRhythmicGenesInHCC.txt",sep = "\t",quote = F,row.names = F,col.names = F)

GainInTumor = na.omit(Ensem_gene[setdiff(rownames(Tumor_rhyGenes),rownames(NAT_rhyGenes)),"SYMBOL"])
write.table(GainInTumor,"~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/GainRhy_commonTCGA/GainRhythmicGenesInHCC.txt",sep = "\t",quote = F,row.names = F,col.names = F)

Hs_Mm <- read.table("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/hs_Mm_hologene.txt")
LostInTumor_mm = Hs_Mm %>% filter(Hs_symbol %in% LostInTumor) %>% select(Mm_symbol,Hs_symbol) %>% unique()


LostRhy_commonTCGA=lapply(Mouse_LostrhyCase,function(x){
  df=na.omit(x)
  data <- intersect(na.omit(x),LostInTumor_mm$Mm_symbol)
  return(data)
})

genelistinConds = list2DF(lapply(LostRhy_commonTCGA, `length<-`, max(lengths(LostRhy_commonTCGA))))
write.xlsx(genelistinConds,"LostRhyGenes_commonHCCLost.xlsx")



bedfile = read.delim("/workspace/rsrch2/common_data/Refgenome/mm10/GenCode.vM25.mm10.Gene_bed/GenCodeVM25.gene.tssup2.5kbDn1kb.padded.filtered.bed",header = F)
names(bedfile)=c("Chr","Start","End","Gene","Score","Strand")


bed_res <- lapply(8:10,function(x){
  genes=na.omit(genelistinConds[,x])
  befile_file = bedfile %>% filter(Gene %in% genes)
  filename = paste0("./LostRhy_Bedfiles_GIGGLE/",names(genelistinConds)[x],".bed")
  #filename=paste0("./GainRhyInCase/Gain_TF.BedFiles/CommonInGT",names(genelist)[x],"condition.bed")
  write.table(befile_file,file = filename,sep = "\t",row.names = F,col.names = F,quote = F)
})


###############Enhancer bed files####
Genebedfile = read.delim("/workspace/rsrch2/common_data/Refgenome/mm10/GenCode.vM25.mm10.Gene_bed/GenCodeVM25.gene.geneBodyUPDn100kb.padded.filtered.bed",header = F)
names(Genebedfile)=c("Chr","Start","End","Gene","Score","Strand")


bed_Enh <- lapply(8:10,function(x){
  genes=na.omit(genelistinConds[,x])
  befile_file = Genebedfile %>% filter(Gene %in% genes)
  filename = paste0("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/LostRhyInCase/LostRhy_commonTCGA/LostRhy_BedfilesMotifs/",names(genelistinConds)[x],".bed")
  # filename = paste0("./GainRhyInCase/Gain_Enhancer.BedFiles/CommonInGT",names(genelist)[x],"condition.bed")
  write.table(befile_file,file = filename,sep = "\t",row.names = F,col.names = F,quote = F)
})



#############################Giggle results summary for lostRhy genes#############
metaInfo = read.delim("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/Cistrom_mm10_tranfac_QC.txt",header = T) %>% select(2,3,5,6)
path="~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/LostRhyInCase/LostRhy_commonTCGA/LostRhy_Bedfiles_GIGGLE/"
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
write.xlsx(df_res,paste0(output,"LostRhyGIGGLE_resforatleat8conditions.xlsx"))


#########################Enhancer motifs results####
homerpath="~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/LostRhyInCase/LostRhy_commonTCGA/LostRhy_BedfilesMotifs/LostRhy_Overlap_Groseq/LostRhy_HomerMotifs_res/"
known_motifsRes <- list.files(path = homerpath,pattern = "knownResults.txt",full.names = T,recursive = T)
Gro_Seq <- known_motifsRes[grep("/Gro",known_motifsRes)]
list_name=gsub(".*/(.*)/.*","\\1",Gro_Seq)

Motif_extract <- function(x){
  data_df <-read.delim(x,header = T,stringsAsFactors = F)
  cond = gsub("GroSeq_","",basename(dirname(x))) 
  data_df$X..of.Target.Sequences.with.Motif = as.numeric(gsub("%","",data_df$X..of.Target.Sequences.with.Motif))
  df_filter <- data_df %>% filter(P.value <= 1e-2 & X..of.Target.Sequences.with.Motif >= 10) %>% 
    mutate(TF =  gsub("(.*)\\(.*\\).*","\\1",gsub("(.*)\\(.*\\).*","\\1",Motif.Name)),
           TF_tyepe = gsub("\\((.*)\\)","\\1",gsub(".*(\\(.*\\)).*\\(.*","\\1",Motif.Name)),
           group=cond
    ) %>% select(Motif.Name,TF,P.value,TF_tyepe,group)
  return(df_filter)
}

TFres_list = lapply(Gro_Seq,Motif_extract)
names(TFres_list)=list_name
TF_res <- do.call(rbind,TFres_list)
TFres_list[[length(TFres_list)+1]]=TF_res
names(TFres_list)[length(TFres_list)]="all_resultsCombined"

write.xlsx(TFres_list,paste0(homerpath,"Groseq_KnownMotifs_res_sheets.xlsx"))




############################################################Gain Rhythmic in HCC###
rm(list = ls())
setwd("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/GainRhy_commonTCGA")

Ensem_gene <- read.delim("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/Ensemble_gene_mappingHsp.txt") %>% 
  group_by(SYMBOL) %>% dplyr::slice(1) %>% ungroup() %>% as.data.frame()
rownames(Ensem_gene)=Ensem_gene$ENSEMBL

Hs_Mm <- read.table("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/hs_Mm_hologene.txt")

TCGA = get(load("~/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/OUT/NATvsTumore.sample.FiT.OUT_Status.RData"))
Liver_NAT=TCGA$`TCGA-LIHC-NAT`$data.fit
Liver_Tumor = TCGA$`TCGA-LIHC-Tumor`$data.fit
Tumor_rhyGenes=Liver_Tumor %>% filter(pval < 0.05) 
NAT_rhyGenes = Liver_NAT %>% filter(pval < 0.05)
GainInTumor = na.omit(Ensem_gene[setdiff(rownames(Tumor_rhyGenes),rownames(NAT_rhyGenes)),"SYMBOL"])
GainInTumor_mm = Hs_Mm %>% filter(Hs_symbol %in% GainInTumor) %>% select(Mm_symbol,Hs_symbol) %>% unique()

Mouse_GainrhyCase <- read.xlsx("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/NonCoreKO.GainRhymic_CaseConditionvsWT_2Sheets.xlsx",sheet = 2) %>%
                     as.list(na.omit(.))


GainRhy_commonTCGA=lapply(Mouse_GainrhyCase,function(x){
  df=na.omit(x)
  data <- intersect(na.omit(x),GainInTumor_mm$Mm_symbol)
  return(data)
})

GaingenelistinConds = list2DF(lapply(GainRhy_commonTCGA, `length<-`, max(lengths(GainRhy_commonTCGA))))
write.xlsx(GaingenelistinConds,"~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/GainRhy_commonTCGA/GainRhyGenes_commonHCCLost.xlsx")

####get bed files for GIGGLE search######
bedfile = read.delim("/workspace/rsrch2/common_data/Refgenome/mm10/GenCode.vM25.mm10.Gene_bed/GenCodeVM25.gene.tssup2.5kbDn1kb.padded.filtered.bed",header = F)
names(bedfile)=c("Chr","Start","End","Gene","Score","Strand")
bed_res <- lapply(8:10,function(x){
  genes=na.omit(GaingenelistinConds[,x])
  befile_file = bedfile %>% filter(Gene %in% genes)
  filename = paste0("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/GainRhy_commonTCGA/GainRhy_Bedfiles_GIGGLE/Gain_",names(GaingenelistinConds)[x],".bed")
  #filename=paste0("./GainRhyInCase/Gain_TF.BedFiles/CommonInGT",names(genelist)[x],"condition.bed")
  write.table(befile_file,file = filename,sep = "\t",row.names = F,col.names = F,quote = F)
})


#####################Enhancer bed files for Motifs
Genebedfile = read.delim("/workspace/rsrch2/common_data/Refgenome/mm10/GenCode.vM25.mm10.Gene_bed/GenCodeVM25.gene.geneBodyUPDn100kb.padded.filtered.bed",header = F)
names(Genebedfile)=c("Chr","Start","End","Gene","Score","Strand")


bed_Enh <- lapply(8:10,function(x){
  genes=na.omit(GaingenelistinConds[,x])
  befile_file = Genebedfile %>% filter(Gene %in% genes)
  filename = paste0("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/GainRhy_commonTCGA/GainRhy_enhancerMotif/Gain_",names(GaingenelistinConds)[x],".bed")
  # filename = paste0("./GainRhyInCase/Gain_Enhancer.BedFiles/CommonInGT",names(genelist)[x],"condition.bed")
  write.table(befile_file,file = filename,sep = "\t",row.names = F,col.names = F,quote = F)
})



########################GIGGLE results for gain Rhythmic genes in HCC and common in Mouse####
metaInfo = read.delim("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/Cistrom_mm10_tranfac_QC.txt",header = T) %>% select(2,3,5,6)
path="~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/GainRhy_commonTCGA/GainRhy_Bedfiles_GIGGLE/"
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
write.xlsx(df_res,paste0(output,"GainRhy_GIGGLE_resforatleat8conditions.xlsx"))



##################Enhancer motifs#####
homerpath="~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/GainRhyInCase/GainRhy_commonTCGA/GainRhy_enhancerMotif/GainRhy_Overlap_GroSeq/GainRhy_HomerMotifs_res/"
known_motifsRes <- list.files(path = homerpath,pattern = "knownResults.txt",full.names = T,recursive = T)
Gro_Seq <- known_motifsRes[grep("/Gro",known_motifsRes)]
list_name=gsub(".*/(.*)/.*","\\1",Gro_Seq)

Motif_extract <- function(x){
  data_df <-read.delim(x,header = T,stringsAsFactors = F)
  cond = gsub("GroSeq_","",basename(dirname(x))) 
  data_df$X..of.Target.Sequences.with.Motif = as.numeric(gsub("%","",data_df$X..of.Target.Sequences.with.Motif))
  df_filter <- data_df %>% filter(P.value <= 1e-2 & X..of.Target.Sequences.with.Motif >= 10) %>% 
    mutate(TF =  gsub("(.*)\\(.*\\).*","\\1",gsub("(.*)\\(.*\\).*","\\1",Motif.Name)),
           TF_tyepe = gsub("\\((.*)\\)","\\1",gsub(".*(\\(.*\\)).*\\(.*","\\1",Motif.Name)),
           group=cond
    ) %>% select(Motif.Name,TF,P.value,TF_tyepe,group)
  return(df_filter)
}

TFres_list = lapply(Gro_Seq,Motif_extract)
names(TFres_list)=gsub("Groseq","",list_name)
TF_res <- do.call(rbind,TFres_list)
TFres_list[[length(TFres_list)+1]]=TF_res
names(TFres_list)[length(TFres_list)]="all_resultsCombined"

write.xlsx(TFres_list,paste0(homerpath,"Groseq_KnownMotifs_res_sheets.xlsx"))



