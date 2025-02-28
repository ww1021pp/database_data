rm(list = ls())
setwd("/workspace/rsrch1/tmp/DataBase_datatable0618/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict")
library(edgeR)
library(data.table)
library(parallel)
library(R.utils)
library(dplyr)

######## to CPM ####

#### Functions #####

#Filter genes with less than 10 counts on average, normalized by library size (TMM) and log2 transform with 0.25 pseudo-count 
run_edgeR = function(k){
  
  posx=which(Proj_ID == k)
  #gtex.sub=as.matrix(NAT_counts[,posx])
  gtex.sub=as.matrix(rawcounts[,posx])
  #rownames(gtex.sub)=paste(gtex$Name,gtex$Description,sep="_")
  gtex.sub=subset(gtex.sub,rowMeans(gtex.sub)>10)
  Ex.y <- DGEList(counts=gtex.sub)
  Ex.y <- calcNormFactors(Ex.y)
  CPM = data.frame(cpm(Ex.y, log=TRUE, prior.count = 0.25))
  names(CPM)=gsub("\\.","-",names(CPM))
  return(CPM)
}

# Remove covariates using a linear regression with age, sex
remove_covariates=function(tt, samp.all, mean.thresh){
  
  tt=subset(tt,rowMeans(tt)>mean.thresh)
  tt=tt[,names(tt)%in%samp.all$SAMPID]
  tt.info=samp.all[names(tt),]
  tt.norm=tt
  
  if(ncol(tt)>20){
    for(j in 1:nrow(tt)){
      for.fit=data.frame(y=as.numeric(tt[j,]), 
                         age=tt.info$AGE, 
                         sex=as.factor(tt.info$SEX))
      
      if(length(unique(for.fit$sex))==1){
        resid=summary(lm(formula = y ~ age  , data = for.fit))$residuals
      }else{
        resid=summary(lm(formula = y ~ age + sex , data = for.fit))$residuals
      }
      
      
      tt.norm[j,]=NA
      tt.norm[j,as.numeric(names(resid))]=resid
    }
    
  }else{
    tt.norm= NULL
  }
  return(tt.norm)
}

spliti= function(x, sp, nb){
  v=sapply(strsplit(x,sp),"[[",nb)
  return(v)
}
#################

dir.create(file.path("./CPM"), showWarnings = FALSE)


#step1 read NAT sample number  gt 25 raw counts###
rawcounts <- read.table("matchedSample_rawCounts.txt",header = T)
colnames(rawcounts)=gsub("\\.","-",colnames(rawcounts))
row_name <- read.table("Gene_name.txt",header = T)
rownames(rawcounts)=row_name$x


sampleInfo <- read.table("MatchSample.metaInfo.txt",header = T) 
rownames(sampleInfo)=sampleInfo$tcga.cgc_sample_id



############convert ensemble ID to Gene Symbol ####
library(org.Hs.eg.db)
gene_ids <-data.frame(
      ENSEMBL=stringr::str_replace(row_name$x, pattern = ".[0-9]+$", replacement = "")
)
rownames(rawcounts)=gene_ids$ENSEMBL
symb <- keys(org.Hs.eg.db, "SYMBOL")
Symbol_ensm<-select(org.Hs.eg.db, symb, "ENSEMBL", "SYMBOL")

tmp = left_join(gene_ids,Symbol_ensm) %>% group_by(ENSEMBL) %>% 
  dplyr::slice(1) %>% ungroup()

rownames(tmp)=tmp$ENSEMBL
write.table(tmp,"Ensemble_gene_mappingHsp.txt",row.names = F,sep = "\t",quote = F)


#############for NAT and Tumor sample ######3
Proj_ID = sampleInfo[colnames(rawcounts),"tcga.gdc_cases.project.project_id"]

uniq_proj <- unique(sampleInfo$tcga.gdc_cases.project.project_id)
Allsample <- sampleInfo %>% 
  mutate(
    SAMPID = tcga.cgc_sample_id,
    PROJID = tcga.gdc_cases.project.project_id,
    CASEID = tcga.cgc_case_id,
    SAMTYPE = tcga.cgc_sample_sample_type,
    SEX = ifelse(tcga.cgc_case_gender == "MALE",1,2),
    AGE = tcga.cgc_case_age_at_diagnosis) %>% 
  dplyr::select(SAMPID,PROJID,CASEID,SAMTYPE,SEX,AGE)
CPM.all= mclapply(uniq_proj, run_edgeR, mc.cores = 16)
names(CPM.all)=uniq_proj
save(CPM.all, file = "./CPM/NATandTumor.CPM_fullSamplegt25.RData")

load("./CPM/NATandTumor.CPM_fullSamplegt25.RData")


CPM.TumorandNAT.norm=mclapply(CPM.all, remove_covariates, Allsample, 1, mc.cores = 16, mc.preschedule = TRUE)


NAT.sample <- sampleInfo %>% filter(tcga.cgc_sample_sample_type == "Solid Tissue Normal") %>%
  mutate(
    SAMPID = tcga.cgc_sample_id,
    PROJID = tcga.gdc_cases.project.project_id,
    CASEID = tcga.cgc_case_id,
    SAMTYPE = tcga.cgc_sample_sample_type,
    SEX = ifelse(tcga.cgc_case_gender == "MALE",1,2),
    AGE = tcga.cgc_case_age_at_diagnosis) %>% 
  dplyr::select(SAMPID,PROJID,CASEID,SAMTYPE,SEX,AGE)



CPM.NAT.Norm=lapply(CPM.TumorandNAT.norm,function(x){
  x=x[,names(x)%in%NAT.sample$SAMPID]
  return(x)
})
save(CPM.NAT.Norm, file = "./CPM/NAT.CPM_fullSamplegt25.RData")

##########Tumor sample########3
Tumor.sample <- sampleInfo %>% filter(tcga.cgc_sample_sample_type != "Solid Tissue Normal") %>%
  mutate(
    SAMPID = tcga.cgc_sample_id,
    PROJID = tcga.gdc_cases.project.project_id,
    CASEID = tcga.cgc_case_id,
    SAMTYPE = tcga.cgc_sample_sample_type,
    SEX = ifelse(tcga.cgc_case_gender == "MALE",1,2),
    AGE = tcga.cgc_case_age_at_diagnosis) %>% 
  dplyr::select(SAMPID,PROJID,CASEID,SAMTYPE,SEX,AGE)


CPM.Tumor.Norm=lapply(CPM.TumorandNAT.norm,function(x){
  x=x[,names(x)%in%Tumor.sample$SAMPID]
  return(x)
})
save(CPM.Tumor.Norm, file = "./CPM/Tumor.CPM_fullSamplegt25.RData")




