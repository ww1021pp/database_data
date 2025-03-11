rm(list=ls())
setwd("~/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/Bed_files/TF_enhancer_drugAble_Res")

####define function###
######
dryplot<-function (dryList, gene1, period = 24) 
{
  gene = Ensem_gene$ENSEMBL[which(Ensem_gene$SYMBOL == gene1)]
  drugcol = allTF$drugCol[which(allTF$Gene == gene1)]

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
      geom_point(size = 0.1, shape = 19) + 
      xlab("Time (h)") + ylab("Log2 normalized counts") + ggtitle(gene1) + 
      scale_x_continuous(breaks = seq(0, period + 6, 6)) + 
      theme_bw(base_size = 8) + theme(aspect.ratio = 1, panel.grid.minor = element_blank(), 
                                      legend.position = "right",
                                      plot.title = element_text(color = drugcol)) + 
      geom_line(aes(x = time, 
                    y = (value), group = group), data = m, position = position_dodge(width = 0.5)) + 
      scale_color_manual(values = c("blue","red"))+
      facet_wrap(~group)
    gg1
  }
}


####load dryR results object###
dryR_res <- readRDS("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/dryR_res/NATsamplegt25.dryR_results.rds")
HCC_dryList = dryR_res$`TCGA-LIHC`
Ensem_gene <- read.delim("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/Ensemble_gene_mappingHsp.txt") %>% group_by(SYMBOL) %>% slice(1) %>% ungroup() %>% as.data.frame()
rownames(Ensem_gene)=Ensem_gene$ENSEMBL

########GRo_seq and TF ChIPseq TF results file###
TFfiles="../ChIPseqGroSeqCombined_TF.xlsx"
TFsheets <- openxlsx::getSheetNames(TFfiles)
TFsheetsList <- lapply(TFsheets,openxlsx::read.xlsx,xlsxFile=TFfiles)
names(TFsheetsList) <- TFsheets
allTF = do.call(rbind,TFsheetsList) %>% select(Gene,col) %>% 
  unique() %>%mutate(drugCol=case_when(
    col == "Aprroved" ~ "red", # target for a approved drug
    col == "Chem" ~ "black", ## not a target for a chemical or drug
    col == "NoChem" ~ "grey" # a target for a drug or Chemical, but not approved 
  )) 
rownames(allTF)=NULL

library(parallel)
for (i in 1:length(TFsheetsList)){
  genes = intersect(toupper(TFsheetsList[[i]]$Gene),toupper(Ensem_gene[rownames(HCC_dryList$parameters),"SYMBOL"]))
  figlist=lapply(genes,function(x){
    p1 <-dryplot(HCC_dryList,gene =x)
  })
  pdf(paste0("./dryR_plots/HCC_TFsPlot_",names(TFsheetsList)[i],".pdf"),height = 2,width = 3)
  for (i in 1:length(figlist)){print(figlist[[i]])}
  dev.off()
}




##################

#####get drug Information###
####based on the plots ##

TFs <- c("CEBPA","NR2F2","BCL6","FOXO1","MYC","FOXO3","EPAS1","ESRRA")

DGinter <- read.delim("~/project/CirDatabase_020824/4DBsets_ChEMBLDGI_Drdugbank_DTC.DTI.txt",header = T) %>% filter(dg_name!="")

DG_TFfilter <- DGinter %>% filter(Gene %in% TFs ) %>% select(Gene,dg_name,dg_group) %>%unique() 


TF_drugFreq=DG_TFfilter %>% mutate(Approved = ifelse(grepl("^[Aa]pproved|_approved",dg_group),1,0)) %>%
  select(Gene,dg_name,Approved) %>% unique() %>%
  group_by(Gene,Approved) %>% summarise(DrugNum=n()) %>% arrange(desc(Approved),desc(DrugNum))
temp = DG_TFfilter %>% select(Gene,dg_name) %>% unique() %>% group_by(Gene) %>% summarise(TotalNum=n()) %>% arrange(TotalNum)
TF_drugFreq$Gene = factor(TF_drugFreq$Gene,levels = temp$Gene)

####get the bar plot of TF with approved drug or chemical ##
TF_drugFreq <- TF_drugFreq %>% mutate(ChemeicalNum = ifelse(Approved == 1,
                                                            log10(DrugNum+1),
                                                            -1*log10(DrugNum+1)))
TF_drugFreq$Approved=factor(TF_drugFreq$Approved)
breaks_values <- pretty(TF_drugFreq$ChemeicalNum)
p2 = ggplot(TF_drugFreq,aes(x = Gene, y = ChemeicalNum, fill = Approved))+
  geom_hline(yintercept = 0)+
  geom_bar(stat = "identity")+
  coord_flip()+
  scale_y_continuous(breaks = breaks_values,
                     labels = abs(breaks_values))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text( color = "black", 
                                    size =8),
        axis.text.y = element_text( color = "black", 
                                    size = 8)
  )+ylab("log10(Chemical Num)")+
  scale_fill_manual(values = c("blue", "red"))
pdf("./LiverTF.inHCC_drugNumer.pdf",width = 3.5, height = 4)
print(p2)
dev.off()



