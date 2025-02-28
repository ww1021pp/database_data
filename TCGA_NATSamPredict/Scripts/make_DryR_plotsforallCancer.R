rm(list = ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/dryR_res/")

library(data.table)
library(dryR)
library(dplyr)
library(parallel)

######
dryplot<-function (dryList, gene, period = 24) 
{
  gene1 = Ensem_gene[gene,"SYMBOL"]
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
                                       legend.position = "right") + 
      geom_line(aes(x = time, 
                                                                                  y = (value), group = group), data = m, position = position_dodge(width = 0.5)) + 
      scale_color_manual(values = c("blue","red"))+
      facet_wrap(~group)
    gg1
  }
}

#########Liver cancer DryR##
dryR_res <- readRDS("NATsamplegt25.dryR_results.rds")

Ensem_gene <- read.delim("../Ensemble_gene_mappingHsp.txt") %>% group_by(SYMBOL) %>% slice(1) %>% ungroup() %>% as.data.frame()
rownames(Ensem_gene)=Ensem_gene$ENSEMBL

library(grid)
library(gridExtra)

for (i in 9:length(dryR_res)){
dryList= dryR_res[[i]]

type=names(dryR_res)[i]
genes = rownames(dryList[["results"]])
gens_pair = Ensem_gene[genes,] %>% group_by(SYMBOL) %>% slice(1)
ifelse(dir.exists(paste0("./plot/",type,"/")),1,dir.create(paste0("./plot/",type,"/")))
figlist=mclapply(gens_pair$ENSEMBL,function(x){
  sym <- Ensem_gene[x,"SYMBOL"]
  p1 <-dryplot(dryList,gene =x)
  ggsave(paste0("./plot/",type,"/",type,"_DryRplot_",sym,".png"),width = 3, height = 3, dpi = 600)
  },mc.cores = 40
)
}






