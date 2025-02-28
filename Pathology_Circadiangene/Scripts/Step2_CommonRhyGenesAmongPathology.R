rm(list=ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/Condition_CirGenesCorrelation/Data_res/Pathology_ChangedGenes/")

PathOlogy_Cir = readRDS("./NineCondition_CircadianList.rds")
NRF_Cir = readRDS("./FiveNRFCircadianList.rds")

PathOlogyList=lapply(PathOlogy_Cir,function(x){
  data=x
  return(data$Gene)
  })

names(PathOlogyList)=gsub("\\.txt","",names(PathOlogyList))
##################Pathology#####

Pathology_binary = t(mtabulate(PathOlogyList))%>% data.frame() %>% rownames_to_column(.,var = "Gene")

GeneinCondition=rowSums(Pathology_binary[,-1])
Pathology_binary$sum=GeneinCondition
Pathology.list = lapply(unique(sort(Pathology_binary$sum)),function(x){
  return(Pathology_binary$Gene[which(Pathology_binary$sum >= x)])
}
)

names(Pathology.list)=paste0("gt",unique(sort(Pathology_binary$sum)))

Pathology_genNumber=cbind(
  V1=names(Pathology.list),
  Freq=as.numeric(sapply(Pathology.list,length))
  ) %>% as.data.frame()

Pathology_genNumber$Gene_number=as.numeric(Pathology_genNumber$Freq)

ggplot(Pathology_genNumber[-1,], aes(x=V1, y=Gene_number,)) + 
  geom_bar(stat = "identity",fill="skyblue")+ 
  geom_text(aes(label = Freq), vjust = 0) + theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2), 
        axis.title = element_text(size = 10),
        plot.title = element_text(size = 10),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(1.5, "pt"))
ggsave("BarplotofPathologySharedCircadianGenes.pdf",height = 3,width = 3)


########################Night Restricted Feeding#######

NRFList=lapply(NRF_Cir,function(x){
  data=x
  return(data$Gene)
})
names(NRFList)=gsub("\\.txt","",names(NRFList))

NRF_binary = t(mtabulate(NRFList))%>% data.frame() %>% rownames_to_column(.,var = "Gene")

NRFGeneinCondition=rowSums(NRF_binary[,-1])
NRF_binary$sum=NRFGeneinCondition
NRF.list = lapply(unique(sort(NRF_binary$sum)),function(x){
  return(NRF_binary$Gene[which(NRF_binary$sum >= x)])
}
)

names(NRF.list)=paste0("NRFgt",unique(sort(NRF_binary$sum)))

NRF_genNumber=cbind(
  V1=names(NRF.list),
  Freq=as.numeric(sapply(NRF.list,length))
) %>% as.data.frame()

NRF_genNumber$Freq=as.numeric(NRF_genNumber$Freq)

ggplot(NRF_genNumber[-1,], aes(x=V1, y=Freq,)) + 
  geom_bar(stat = "identity",fill="skyblue")+ 
  geom_text(aes(label = Freq), vjust = 0) + theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2), 
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(1.5, "pt"))
ggsave("BarplotofNRFSharedCircadianGenes.pdf",height = 3,width = 3)


SharedGenes = list(
  PathOlogy = list2DF(lapply(Pathology.list, `length<-`, max(lengths(Pathology.list)))),
  NRF = list2DF(lapply(NRF.list, `length<-`, max(lengths(NRF.list))))
  )

write.xlsx(SharedGenes,"CommonCirGenesinPathologyOrNRF_2Sheets.xlsx")


library(ggVennDiagram)
# install.packages("ggplot2")
library(ggplot2)

figs=lapply(2:7, function(x){
  na = names(Pathology.list[[x]])
  plot_list = append(NRF.list[-1],Pathology.list[x])
  ggVennDiagram(plot_list, color = "black", lwd = 0.8, lty = 1,label = "count",
                label_alpha=0,label_size = 4) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF")+
    theme(legend.position = "none") 
})


plot_list <- lapply(figs, function(p) p + theme(plot.margin = margin(5, 5, 5, 5)))

pdf("my_plots.pdf", width = 11, height = 10)
marrangeGrob(grobs = plot_list, ncol = 2, nrow = 3)
dev.off()
