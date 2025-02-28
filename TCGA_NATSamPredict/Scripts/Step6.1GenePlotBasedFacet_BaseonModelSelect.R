rm(list = ls())
setwd("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict")
library(vroom)
library(ggplot2)
library(tibble)
library(gridExtra)
library(ggrepel)
library(dplyr)
########

###### plot file based through facet#####
Plot_profiles<-function(OUT,tissues,geni, plt.hist=FALSE, R.plot=FALSE, val="pval", text_size=10, title_add="", always_plot_fit=F, vcut=0.05, xl=NULL, period=24){
  p2=period/2
  nmz=names(OUT)

  tissuex= nmz[grepl(tissues,nmz)]
  pltlist=list()
  genepltList=list()
    for(i in tissuex){ 
    group=gsub(paste0(tissues,"-"),"",i)
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
    gene = Ensem_gene$ENSEMBL[which(Ensem_gene$SYMBOL == geni)]
    geni_com = intersect(gene,gene.list)
      fg=geni_com
      geneSymbol=names(geni)[which(geni == fg)]
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
  labels_a=unique(paste0(geni, '\npval ', signif(geneplt$pval[],2)))
  geneplt$col=factor(paste0(geneplt$type,"_",geneplt$group))
  geneplt$group=factor(geneplt$group, levels = c("NAT","Tumor"), 
                       labels = labels_a )
  
        pltlist=ggplot()+geom_point(data=geneplt[which(geneplt$type=="data"),],aes(x=time, y=gene, color=col, size=type,alpha=type))+
          geom_line(data=geneplt[which(geneplt$type=="fit"),],aes(x=time,y=gene,color=col),linewidth=0.5)+
          theme_bw()+ xlim(0,24)+
          labs(x=xl ,y="Centred expression [log2]")+
          theme(legend.position = "none",text=element_text(size=text_size))+scale_size_manual(values=c(0,0.5,0))+
          scale_alpha_manual(values=c(1,0.8,1))+scale_color_manual(values=c("#008ECC","#008ECC","blue","red"))+facet_wrap(~group,scales="free")
     
    return(pltlist)
    } 


load("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/OUT/NATvsTumore.sample.FiT.OUT_Status.RData")

Ensem_gene <- read.delim("Ensemble_gene_mappingHsp.txt") %>% group_by(SYMBOL) %>% dplyr::slice(1) %>% ungroup() %>% as.data.frame()
rownames(Ensem_gene)=Ensem_gene$ENSEMBL

rownames(OUT.NT$`TCGA-LUSC-NAT`$data.fit)[1:10]

load("./OUT/SS_SampleTypeModelselect.RData")
Liver = SS$`TCGA-LIHC`
tmp=Liver[fg,ncol(Liver)-20,ncol(Liver)]

#gene_inf <- read.delim("/workspace/rsrch2/panpanliu/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/Scripts/CRG_EN.txt",header = T,sep=" ")
#geni = gene_inf$x;names(geni)=rownames(gene_inf)
#geni[13]="ENSG00000173153"
#names(geni)[8]="BMAL1"
#names(geni)[13]="ESRRA"

tissues = unique(gsub("(.*-.*)-.*","\\1",names(OUT.NT)))

for (i in tissues){
  genes = rownames(OUT.NT[[paste0(i,"-NAT")]]$data.fit)
  gens_pair = Ensem_gene[genes,] %>% group_by(SYMBOL) %>% dplyr::slice(1)
  gens_pair=gens_pair[1:100,]
  ifelse(dir.exists(paste0("./plot/",i,"/")),1,dir.create(paste0("./plot/",i,"/")))
  figlist=mclapply(gens_pair$ENSEMBL,function(x){
  sym <- Ensem_gene[x,"SYMBOL"]
  p1 <-Plot_profiles(geni = sym,OUT = OUT.NT, tissues = i)
  ggsave(paste0("./plot/",i,"/",i,"_ModelSelect_",sym,".png"),width = 4, height = 2, dpi = 600)
},mc.cores = 40
)
}


