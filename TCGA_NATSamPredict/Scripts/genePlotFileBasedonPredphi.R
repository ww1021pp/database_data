rm(list = ls())
setwd("/workspace/rsrch1/tmp/DataBase_datatable0618/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict")
library(vroom)
library(ggplot2)
library(tibble)
library(gridExtra)
library(ggrepel)
########

Plot_profiles<-function(OUT, tissues, geni, plt.hist=FALSE, R.plot=FALSE, val="pval", text_size=10, title_add="", always_plot_fit=F, vcut=0.2, xl=NULL, period=24){
  p2=period/2
  nmz=names(OUT)
 
  tissuex= intersect(nmz, tissues)
  if(length(tissuex)==0){
    tissuex=NULL
    for (nm in tissues) {
      ttix=names(OUT)[startsWith(names(OUT),nm)]
      tissuex=c(tissuex,ttix)
    }
  }
  for(i in tissuex){ 
    pltlist=list()
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
    geni_com = intersect(geni,gene.list)
    for(k in 1:length(geni_com)){
      fg=geni_com[k]
      geneSymbol=names(geni)[which(geni == fg)]
      idx=match(fg, gsub("\\|.*$","",gsub("^.*_", "",rownames(out$E))))
      if(is.na(idx)){
        cat("gene ", fg, " not present in ", i, "\n")
        stop()
      }
      ts=(1:1000)*pi/500
      pred=GE[fg]*Conj(complex(argument = (ts*24/period)))#+mus[fg]
      fit=tibble(gene=Re(pred), phase=ts,type="fit", time=ts*12/pi, sd=0, w=0)#, conf=1)
      geneplt=tibble(gene=as.matrix(out$E)[idx,],phase=out$phi,type="data", time=out$phi*12/pi, sd=0, w=0)#, conf=out$score)
      aver=sdev=1:12
      phibinz=c(0:11)*pi/6
      for (bin in 1:(length(phibinz))) {
        phizb=which(out$phi>phibinz[bin])
        phizs=which(out$phi<phibinz[bin]+pi/6)
        phiz=intersect(phizb, phizs)
        aver[bin]=mean(as.matrix(out$E)[idx,phiz], na.rm = T)
        sdev[bin]=sd(as.matrix(out$E)[idx,phiz], na.rm=T)/sqrt(length(phiz))
      }
      tibaver=tibble(gene=aver, phase=(phibinz+pi/12), type="avg", time=(phibinz+pi/24)*12/pi, sd=sdev, w= 0.2)
      if (always_plot_fit) geneplt=rbind(geneplt, tibaver, fit)
      else{
        if(fut[[val]][idx]<vcut) geneplt=rbind(geneplt, tibaver, fit)
        else geneplt=rbind(geneplt, tibaver)
      }
      if(! (val%in%colnames(fut))){val="pval"}
      if(k==1) {
        phi.tb=tibble(phi=inf.phi*12/pi)
        brk=seq(0, 24, 1)
        
        if(!is.character(title_add)) ttl=NULL
        else ttl=paste(i, title_add)
        
        if(plt.hist) {pltlist[[l]]=ggplot(phi.tb, aes(x=phi))+geom_histogram(breaks=brk,fill="#008ECC")+
          theme_minimal()+ lims(x=c(0,24))+labs(x=xl,y="Frequency", title=ttl) +
          theme(text=element_text(size=text_size))
        l=l+1
        pltlist[[l]]=ggplot(geneplt)+geom_point(aes(x=time, y=gene, color=type, size=type,alpha=type))+
          labs(x=xl,y="Centred expression [log2]", subtitle=paste(geneSymbol, val, signif(fut[[val]][idx],2)))+theme_minimal()+ lims(x=c(0,24))+#c(-5,5))+
          theme(legend.position = "none",text=element_text(size=text_size))+scale_size_manual(values=c(3,1,0.00007))+scale_alpha_manual(values=c(1,0.8,1))+scale_color_manual(values=c("#008ECC","#B0DFE5","#111E6C","#0B49BD"))
        #+geom_errorbar(aes(x=time, y=gene, ymin=gene-sd, ymax=gene+sd,colour="errorbar", width=w))
        }
        else{ pltlist[[l]]=ggplot(geneplt)+geom_point(aes(x=time, y=gene, color=type, size=type,alpha=type))+
          labs(x=xl,y="Centred expression [log2]", title=ttl, subtitle=paste(geneSymbol, val, signif(fut[[val]][idx],2)))+
          theme_minimal()+ lims(x=c(0,24))+#c(-5,5))+
          theme(legend.position = "none",text=element_text(size=text_size))+scale_size_manual(values=c(3,1,0.00007))+scale_alpha_manual(values=c(1,0.8,1))+scale_color_manual(values=c("#008ECC","#B0DFE5","#111E6C","#0B49BD"))
        #+geom_errorbar(aes(x=time, y=gene, ymin=gene-sd, ymax=gene+sd,colour="errorbar", width=w))
        }
        l=l+1
      }
      else {pltlist[[l]]=ggplot(geneplt)+geom_point(aes(x=time, y=gene, color=type, size=type))+
        labs(x=xl ,y="Centred expression [log2]", subtitle=paste(geneSymbol, val, signif(fut[[val]][idx], 2)))+theme_minimal()+ lims(x=c(0,24))+#c(-5,5))+
        theme(legend.position = "none", text=element_text(size=text_size))+scale_size_manual(values=c(3, 1, 0.00007))+scale_color_manual(values=c("#008ECC","#B0DFE5","#111E6C","#0B49BD"))
      #+geom_errorbar(aes(x=time, y=gene, ymin=gene-sd, ymax=gene+sd, color="errorbar",width=w))
      l=l+1}
    }
    if(length(pltlist)==1){f.plot=pltlist[[1]]}
    else f.plot=do.call(grid.arrange,c(grobs=pltlist, as.table=FALSE,ncol = 3 ))
    print(f.plot)
    if(R.plot) return(f.plot)
    }
  
}


load("./OUT/OUT_ALL.NAT.RData")
gene_inf <- read.delim("../Scripts/CRG_EN.txt",header = T,sep=" ")
geni = gene_inf$x;names(geni)=rownames(gene_inf)

pdf("./plot/CoreClockRelatedGene.pdf")
tissuex=c("TCGA-LIHC", "TCGA-BRCA")
geni = geni
Plot_profiles(OUT, tissuex, geni, val="qval")
dev.off()


