rm(list = ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/dryR_res/")

library(data.table)
library(dryR)
library(dplyr)

tp <- read.delim("../OUT/MatchedNAT_sampleTimePoints.txt",header = T) %>% arrange(phi)
rownames(tp)=tp$CaseID

load("../CPM/NAT.CPM_fullSamplegt25.RData")
load("../CPM/Tumor.CPM_fullSamplegt25.RData")

#####define function ####
OUT_CPM = function(x,CPM){
  res=list()
  NAT = CPM[[x]]
  tmp=colnames(NAT)
  names(tmp)=gsub("(.*-.*-.*)-.*","\\1",names(NAT))
  res$timepoints = tp[intersect(tp$CaseID,names(tmp)),"phi"]
  tmp1 = tmp[intersect(tp$CaseID,names(tmp))]
  res$Expr = NAT[,tmp1]
  return(res)
} 

###################### define dryR function for each project ID ####
DryR_NATvsTumor = function(x){
  data = cbind(NAT_CPMtp[[x]]$Expr,Tumor_CPMtp[[x]]$Expr)
  group=rep(c("NAT","Tumor"),each=ncol(NAT_CPMtp[[x]]$Expr))
  time=c(NAT_CPMtp[[x]]$timepoints,Tumor_CPMtp[[x]]$timepoints)
  dryList = drylm(data,group,time)
  return(dryList)
  }

#######################dry plot function with gene name####

dryplot<-function (dryList, gene1, period = 24) 
{
  gene = Ensem_gene$ENSEMBL[which(Ensem_gene$SYMBOL == gene1)]
  if(gene %in% rownames(dryList[["results"]])){
  normal = FALSE
  if ("ncounts" %in% names(dryList)) {
    vsd = log2(dryList[["ncounts"]] + 1)
  }
  if ("values" %in% names(dryList)) {
    vsd = dryList[["values"]]
    normal = T
  }
  parameters = dryList[["parameters"]][, grep("^mean|^a_|^b_|^amp|^phase|^relamp", 
                                              colnames(dryList[["parameters"]]))]
  ID = rownames(dryList[["results"]])[grep(paste0("^", gene, 
                                                  "$"), rownames(dryList[["results"]]))]
 
  d = vsd[ID, ]
  d = reshape2::melt(d)
  d$group = dryList[["group"]]
  d$time = as.numeric(dryList[["time"]])
  d$time = d$time%%period
  suppressWarnings({
    d <- Rmisc::summarySE(d, measurevar = "value", groupvars = c("time", 
                                                                 "group"))
  })
  v = seq(0, period, round(period/24, 0))
  fit_d_0 = parameters[which(rownames(parameters) == ID), grep("mean", 
                                                               colnames(parameters))]
  fit_d_1 = parameters[which(rownames(parameters) == ID), grep("a_", 
                                                               colnames(parameters))]
  fit_d_2 = parameters[which(rownames(parameters) == ID), grep("^b_", 
                                                               colnames(parameters))]
  fit_d_0[is.na(fit_d_0)] = 0
  fit_d_1[is.na(fit_d_1)] = 0
  fit_d_2[is.na(fit_d_2)] = 0
  m = data.frame(v)
  dd = data.frame(v)
  dd$v = v
  fit_values = function(x, n) {
    as.numeric((fit_d_0[n] + fit_d_1[n] * cos(2 * pi * x/period) + 
                  fit_d_2[n] * sin(2 * pi * x/period)))
  }
  for (u in 1:length(unique(d$group))) {
    m[, u + 1] = NA
    m[, u + 1] = apply(dd, 1, fit_values, u)
  }
  m = m[, -1]
  colnames(m) = unique(dryList[["group"]])
  m = reshape2::melt(m, , id.vars = NULL)
  m$time = rep(v, length(unique(d$group)))
  colnames(m) = c("group", "value", "time")
  if (normal == FALSE) {
    m$value[which(m$value < 0)] = 0
  }
  gg1 = ggplot(d, aes(x = time, y = value, group = group, color = group)) + 
    #geom_errorbar(aes(ymin = value - se, ymax = value + se), 
  #                width = 0.4) + 
     geom_point(size = 0.5, shape = 19) + 
    xlab("Time (h)") + ylab("Log2 normalized counts") + ggtitle(gene1) + 
    scale_x_continuous(breaks = seq(0, period + 6, 6)) + 
    theme_bw(base_size = 10) + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), 
                                     legend.position = "right") + 
    geom_line(aes(x = time, y = (value), group = group), data = m, position = position_dodge(width = 0.5)) + 
    scale_color_manual(values = c("blue","red"))+
    facet_wrap(~group)
  #gg1
  ggsave(paste0(output_dir, "plot.pdf"), plot = my_plot, width = 6, height = 4)
}
}


NAT_CPMtp <- lapply(names(CPM.NAT.Norm),OUT_CPM,CPM.NAT.Norm)
names(NAT_CPMtp)=names(CPM.NAT.Norm)

Tumor_CPMtp <- lapply(names(CPM.Tumor.Norm),OUT_CPM,CPM.Tumor.Norm)
names(Tumor_CPMtp)=names(CPM.NAT.Norm)


dryR_res = lapply(1:length(names(CPM.Tumor.Norm)),DryR_NATvsTumor)
names(dryR_res)=names(CPM.NAT.Norm)
saveRDS(dryR_res,"NATsamplegt25.dryR_results.rds")

#########Liver cancer DryR##
dryR_res <- readRDS("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/dryR_res/NATsamplegt25.dryR_results.rds")
HCC_dryList = dryR_res$`TCGA-LIHC`
temp = HCC_dryList$parameters %>% mutate(gene=rownames(.)) %>% 
  group_by(chosen_model) %>% sample_n(5) %>% ungroup()

Ensem_gene <- read.delim("../Ensemble_gene_mappingHsp.txt") %>% group_by(SYMBOL) %>% slice(1) %>% ungroup() %>% as.data.frame()
rownames(Ensem_gene)=Ensem_gene$ENSEMBL
geneset = Ensem_gene[temp$gene,"SYMBOL"]

CRG_gene <- read.table("../../Scripts/CRG_EN.txt")
genes = rownames(CRG_gene)
genes[8]="BMAL1"

figlist=list()
figlist=lapply(genes,function(x){
    p1 <-dryplot(HCC_dryList,gene1 =x)
    return(p1)}
)

library(grid)
library(gridExtra)

pdf("./plot/CoreClockRelatedGene_dryRPlots.pdf",height = 3,width = 4)
for (i in 1:length(figlist)){print(figlist[[i]])}
dev.off()

ggsave(filename = plot.path, plot = Plot(df2, gene.tmp,p_value), width = 3, height = 3, dpi = 600)







