Here is a scripts to plot the genes expression in TCGA tumor and NAT samples

######
```
rm(list = ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict")
library(vroom)
library(ggplot2)
library(tibble)
library(gridExtra)
library(ggrepel)
library(dplyr)
library(parallel)
########

###### plot file based through facet#####
Plot_profiles<-function(OUT,tis.tmp,geni.tmp, plt.hist=FALSE, R.plot=FALSE, val="pval", text_size=10, title_add="", always_plot_fit=F, vcut=0.05, xl=NULL, period=24){
  p2=period/2
  nmz=names(OUT)
  
  tissuex= nmz[grepl(tis.tmp,nmz)]
  pltlist=list()
  genepltList=list()
  for(i in tissuex){ 
    group=gsub(paste0(tis.tmp,"-"),"",i)
    l=1
    sz=7
    ang=90
    hjst=1
    out=OUT[[i]]
    fit=out$data.fit
    fut=fit
    inf.phi=out$phi
    exprx=out$E
    gene.list=gsub("\\|.*$","",gsub("^.*_", "",fit$genes))
    full=fit[,c("amp","a","b","R2", "genes")]
    full$genz=gsub("\\|.*$","",gsub("^.*_", "",full$genes))
    full$R=sqrt(full$a^2+full$b^2)
    GE=fit[,c("a","b", "genes")]
    GL=GE$genes
    GE=complex(real = GE[,"a"], imaginary =  GE[,"b"])
    names(GE)=gsub("\\|.*$","",gsub("^.*_", "",GL))
    if(is.null(xl)) xl="CHIRAL ZT [h]"
    gene = Ensem_gene$ENSEMBL[which(Ensem_gene$SYMBOL == geni.tmp)]
    geni_com = intersect(gene,gene.list)
    fg=geni_com
    geneSymbol=geni.tmp
    idx=match(fg, gsub("\\|.*$","",gsub("^.*_", "",rownames(out$E))))
    if(is.na(idx)){
      cat("gene ", fg, " not present in ", i, "\n")
      stop()
    }
    ts=(1:1000)*pi/500
    pred=GE[fg]*Conj(complex(argument = (ts*24/period)))#+mus[fg]
    fit=tibble(gene=Re(pred), phase=ts,type="fit", time=ts*12/pi, sd=0, w=0,group=group)#, conf=1)
    geneplt=tibble(gene=as.matrix(out$E)[idx,],phase=out$phi,type="data", time=out$phi*12/pi, sd=0, w=0,group=group)#, conf=out$score)
    aver=sdev=1:12
    phibinz=c(0:11)*pi/6
    for (bin in 1:(length(phibinz))) {
      phizb=which(out$phi>phibinz[bin])
      phizs=which(out$phi<phibinz[bin]+pi/6)
      phiz=intersect(phizb, phizs)
      aver[bin]=mean(as.matrix(out$E)[idx,phiz], na.rm = T)
      sdev[bin]=sd(as.matrix(out$E)[idx,phiz], na.rm=T)/sqrt(length(phiz))
    }
    tibaver=tibble(gene=aver, phase=(phibinz+pi/12), type="avg", time=(phibinz+pi/24)*12/pi, sd=sdev, w= 0.2,group=group)
    if(fut[[val]][idx] <= vcut) geneplt=rbind(geneplt, tibaver, fit)
    else {
      fit=tibble(gene=mean(as.matrix(out$E)[idx,]), phase=ts,type="fit", time=ts*12/pi, sd=0, w=0,group=group)#, conf=1)
    }
    geneplt=rbind(geneplt, tibaver,fit) %>% mutate(pval=fut[[val]][idx])
    genepltList[[i]]=rbind(genepltList[[i]],geneplt)
  }
  geneplt=do.call(rbind,genepltList)
  labels_a=unique(paste0(geni.tmp, '\npval ', signif(geneplt$pval[],2)))
  geneplt$col=factor(paste0(geneplt$type,"_",geneplt$group))
  geneplt$group=factor(geneplt$group)
  
  Tumor_data = geneplt %>% filter(group == "Tumor")
  Tumor_pvalue=format(Tumor_data$pval[1], scientific = TRUE,digits = 2)
  NAT_data = geneplt %>% filter(group=="NAT")
  NAT_pvalue=format(NAT_data$pval[1], scientific = TRUE,digits = 2)
  
  theme <- theme_bw() +
    theme(#axis.line = element_blank(),  
          axis.text = element_text(color = "black", size = 10),
          title = element_text(size=12, colour = "black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.8), 
          axis.title = element_text(size = 10,color = "black"),
          plot.title = element_text(size = 10,color = "black"))
  
  
  eq <- substitute(
    expr =
      paste(
        italic("p"),
        " = ",
        NAT_pvalue,
      ),
    env = base::list(NAT_pvalue=NAT_pvalue))
  
  plt_1 <- ggplot()+
    geom_point(data=NAT_data[which(NAT_data$type=="data"),],aes(x=time, y=gene, alpha=type),color="grey30",size=0.5)+
    geom_line(data=NAT_data[which(NAT_data$type=="fit"),],aes(x=time,y=gene),linewidth=0.5,color="blue")+
    scale_x_continuous(limits = c(0,24),breaks=seq(0,24,6))+expand_limits(x=0)+
    labs(x="" ,y="",title = "NAT",subtitle = eq)+theme+
    theme(legend.position = "none") +scale_alpha_manual(values = c(1, 0.8))

  eq_Tu <- substitute(
    expr =
      paste(
        italic("p"),
        " = ",
        Tumor_pvalue,
      ),
    env = base::list(Tumor_pvalue=Tumor_pvalue))
  
  plt_2 <- ggplot()+
    geom_point(data=Tumor_data[which(Tumor_data$type=="data"),],aes(x=time, y=gene, alpha=type),color="grey30",size=0.5)+
    geom_line(data=Tumor_data[which(Tumor_data$type=="fit"),],aes(x=time,y=gene),linewidth=0.5,color="red")+
    theme_bw()+ scale_x_continuous(limits = c(0,24),breaks=seq(0,24,6))+expand_limits(x=0)+
    labs(x="" ,y="",title = gsub("TCGA-","",tis.tmp),subtitle = eq_Tu)+theme+
    theme(legend.position = "none")+scale_alpha_manual(values = c(1, 0.8))
  
  svg(paste0(plotpath,"Transcripts_",tis.tmp,"_",geni.tmp,".svg"), width = 4, height = 3) # Open a new pdf file
   grid.arrange(plt_1, plt_2, nrow = 1,
               bottom = textGrob("Time of day(h)", gp = gpar(fontsize = 12),vjust=-1),
               left = textGrob("log2 normalized expression", rot = 90, gp = gpar(fontsize = 12),vjust=1.5),
               top=textGrob(geni.tmp,  gp = gpar(fontsize = 12,fontface = "bold"),just = "center")
               )
  dev.off()
 
} 



load("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/OUT/NATvsTumore.sample.FiT.OUT_Status.RData")

Ensem_gene <- read.delim("Ensemble_gene_mappingHsp.txt") %>% group_by(SYMBOL) %>% dplyr::slice(1) %>% ungroup() %>% as.data.frame()
rownames(Ensem_gene)=Ensem_gene$ENSEMBL
plotpath="~/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/plot/Figures_forDatabase/"
rownames(OUT.NT$`TCGA-LUSC-NAT`$data.fit)[1:10]

gene_inf <- read.delim("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/Scripts/CRG_EN.txt",header = T,sep=" ")
geni = gene_inf$x;names(geni)=rownames(gene_inf)
geni[13]="ENSG00000173153"
names(geni)[8]="BMAL1"
names(geni)[13]="ESRRA"
tissues = unique(gsub("(.*-.*)-.*","\\1",names(OUT.NT)))
for (i in tissues){
  genes = rownames(OUT.NT[[paste0(i,"-NAT")]]$data.fit)
  gens_pair = Ensem_gene[genes,] %>% filter(SYMBOL %in% names(geni)) %>% group_by(SYMBOL) %>% dplyr::slice(1)
  #gens_pair=gens_pair[1:100,]
  #ifelse(dir.exists(paste0("./plot/",i,"/")),1,dir.create(paste0("./plot/",i,"/")))
  figlist=mclapply(gens_pair$ENSEMBL,function(x){
    sym <- Ensem_gene[x,"SYMBOL"]
    p1 <-Plot_profiles(geni.tmp = sym,OUT = OUT.NT, tis.tmp = i)
   },mc.cores = 10
  )
}

```
![image](https://github.com/user-attachments/assets/9a8c29d5-bc9b-40b9-9796-bb6c9529ac83)


