rm(list=ls())
setwd("/workspace/rsrch1/tmp/DataBase_datatable0618/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict")
#devtools::install_github("naef-lab/CHIRAL/Pkg/CHIRAL")
require(CHIRAL)

head(muscle_exon)      
dim(muscle_exon)
library(lmtest)
library(parallel)
library(R.utils)
library(dplyr)

###################define function#######

CPM_to_E<-function(CPM.all.norm, min.samp=20, sep=NULL, samp=NULL){
  E.matrix=list()
   if(!is.null(samp)){
    if(tolower(sep)=="sex"){
      male=samp$sub.id[samp$SEX==1]
      female=samp$sub.id[samp$SEX==2]
      for(name in names(CPM.all.norm)){
        E=list()
        E$E=CPM.all.norm[[name]]
        E$E=E$E[,!is.na(E$E[1,])]
        if(!is.null(E$E)){
          E$tissue=name
          ml=fm=list()
          colnames(E$E)=spliti(colnames(E$E), "\\.", 2)
          ml$E=E$E[,colnames(E$E) %in% male]
          ml$type="GTEX-male"
          ml$tissue=name
          if(ncol(ml$E)>=min.samp){E.matrix[[paste(name, "Male", sep="-")]]=ml}
          fm$E=E$E[,colnames(E$E) %in% female]
          fm$type="GTEX-female"
          fm$tissue=name
          if(ncol(fm$E)>=min.samp){E.matrix[[paste(name, "Female", sep="-")]]=fm}
        }
      }
    }
    if( tolower(sep)=="age"){
      young=samp$sub.id[samp$AGE>60]
      old=samp$sub.id[samp$AGE<50]
      for(name in names(CPM.all.norm)){
        E=list()
        E$E=CPM.all.norm[[name]]
        E$E=E$E[,!is.na(E$E[1,])]
        if(!is.null(E$E)){
          E$tissue=name
          ml=fm=list()
          ml$E=E$E[,colnames(E$E) %in% young]
          ml$type="GTEX-old"
          ml$tissue=name
          if(ncol(ml$E)>=min.samp){E.matrix[[paste(name, "Young", sep="-")]]=ml}
          fm$E=E$E[,colnames(E$E) %in% old]
          fm$type="GTEX-young"
          fm$tissue=name
          if(ncol(fm$E)>=min.samp){E.matrix[[paste(name, "Old", sep="-")]]=fm}
        }
      }
    }
  }
  if(length(names(E.matrix))==0){
    for(name in names(CPM.all.norm)){
      E=list()
      E$E=CPM.all.norm[[name]]
      E$E=E$E[,!is.na(E$E[1,])]
      if(!is.null(E$E)){
        if(ncol(E$E)>min.samp){
          E$projID=name
          gene="ENSG00000049246"
          E.matrix[[name]]=E}
      }
    }
  }

  return(E.matrix)
}
  
infer_l=function(k,clockgenes=NULL){
  tix=k$projID
  v=k$E
  out=CHIRAL(v, 500, clockgenes = clockgenes,standardize = T) 
  out$tissue=tix
  return(out)
}

#Harmonic regression
harm_reg<-function(x, t, period){
  n=length(x)
  fit0=lm(x~1)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit1=lm(x~c+s)
  mu=coef(fit1)[1]
  a=coef(fit1)[2]
  b=coef(fit1)[3]
  p.val=lrtest(fit1, fit0)$Pr[2]
  amp=2*sqrt(a^2+b^2)
  phase=atan2(b,a)%%(2*pi)
  phase=period*phase/(2*pi)
  return(c(pval=p.val,phase=phase,amp=amp,mu=mu,a=a,b=b, period=period, R2=summary(fit1)$r.squared))
}

#Fit each gene with Harmonic regression in all tissues
Fit_OUT<-function(OUT,period=24, NA5=T, N.cores){
  for(i in names(OUT)){ 
    E=OUT[[i]]$E
    de=sweep(E,1,rowMeans(E),FUN="-")
    if(NA5) E[abs(de)>5]=NA
    phase=OUT[[i]]$phi
    if(length(phase)<24){
      OUT[[i]]=NULL
      next
    }
    
    genes=as.list(rownames(E))
    names(genes)=genes
    dat.fit=mclapply(genes,function(x) harm_reg(as.numeric(E[x,]), 12*as.numeric(phase)/pi, period=period),mc.cores = N.cores)
    dat.fit=do.call(rbind,dat.fit)
    dat.fit=as.data.frame(dat.fit)
    
    dat.fit$qval=p.adjust(dat.fit$pval, "BH")
    OUT[[i]]$data.fit=data.frame(dat.fit,E,genes=rownames(E))
    colnames(OUT[[i]]$data.fit)=c("pval","phase", "amp","mu", "a","b", "period", "R2","qval", colnames(E), "genes")
    OUT[[i]]$E=OUT[[i]]$data.fit[,colnames(E)]
  }
  return(OUT)
}


#Break inveriances intrinic of the method
Set_OUT<-function(OUT){
  for(i in names(OUT)){
    out=OUT[[i]]
    out=order.out.setgene.hc(out)
    sampz=gsub("\\..*$","", gsub("GTEX.","", colnames(out$E)))
    colnames(out$E)=sampz
    names(out$phi)=sampz
    OUT[[i]]=out
  }
  return(OUT)
}
#Break the rotational inveriance
order.out.setgene.hc<-function(out,gene="ENSG00000049246"){
  out$has.been.flipped=0
  d.fit=out[["data.fit"]]
  A=complex(real=d.fit$a, imaginary = d.fit$b)
  names(A)=gsub("\\|.*$","",gsub("^.*_", "",d.fit$genes))
  AA=order.from.hc.ref(A)
  if(sum(Im(A)*Im(AA))<0){
    out[["phi"]]=(-out[["phi"]])%%(2*pi)
    out$has.been.flipped=1
  }
  nz=complex(argument=(pi*(1-1/4)-Arg(AA["ENSG00000049246"])))
  nr=rep(nz,length(AA))
  AA=AA*nr
  out[["phi"]]=(out[["phi"]]+Arg(nz))%%(2*pi)
  d.fit$b=Im(AA)
  d.fit$a=Re(AA)
  d.fit$phase=(Arg(AA)%%(2*pi))*12/pi
  out[["data.fit"]]=d.fit
  return(out)
}

#Break the time direction invariance using human muscle refence
order.from.hc.ref<-function(x){
  sx=x
  #hard coded reference
  full.ref=c(-0.43043911+0.02646768i, -0.21520829-0.08481562i, -0.18976262-0.06323798i, -0.26402067+0.03233651i, 
             -0.15460540-0.17813898i,  0.34503079+0.12988785i,  0.39656252+0.03842312i, -0.30343360+0.29074942i,  0.06226694-0.12542420i, -0.06159846-0.04409744i, -0.18018038-0.01939843i)
  names(full.ref)=c("ENSG00000105516"   ,  "ENSG00000049246"  ,  "ENSG00000167074"   ,  "ENSG00000174738"  ,   "ENSG00000132326"  ,  "ENSG00000170485" ,  "ENSG00000133794" ,  "ENSG00000126368"  , "ENSG00000008405"  ,  "ENSG00000121671",  "ENSG00000179094")
  
  shifts=c(1:1000)*pi/500
  ts=complex(real=cos(shifts), imaginary=sin(shifts))
  common=intersect(names(full.ref), names(x))
  ref.mat=full.ref[common]%o%ts
  x=x[common]
  x=x/(sqrt(sum(Re(x*Conj(x)))))
  gen.p.scal=max(Re(t(ref.mat)%*%Conj(x))) #remember the inversion 
  inv.p.scal=max(Re(t(ref.mat)%*%x))
  
  if(gen.p.scal>inv.p.scal){return(sx)}
  return(Conj(sx))
}

#########get the Donors information##
Unique_donors<-function(OUT){
  all_people=NULL
  for(i in names(OUT)){
    sampz=colnames(OUT[[i]]$E)
    all_people=c(all_people, sampz)
  }
  people=unique(all_people)
  return(people)
}


#Set all the TIPs in a matrix
Create_phi_matrix<-function(OUT, people){
  tix.study=names(OUT)
  phi_mat=matrix(nrow = length(people), ncol= length(names(OUT)))
  dimnames(phi_mat)=list(people, names(OUT))
  for(i in names(OUT)){
    out=OUT[[i]]
    sampz=colnames(out$E)
    phi_mat[sampz,i]=out$phi
  }
  phi_study=phi_mat[,tix.study]
  return(phi_study)
}


###############main function########

####################For TCGA NAT sample predict###
load("./CPM/NAT.CPM_fullSamplegt25.RData")
CPM.all.norm=CPM.NAT.Norm
out_list<-mclapply(CPM.all.norm,CHIRAL, clockgenes = CRG_ens)  

gene_inf <- read.delim("../Scripts/CRG_EN.txt",header = T,sep=" ")
E=CPM_to_E(CPM.all.norm)

OUT=mclapply(E, infer_l, gene_inf$x, mc.cores=16)

OUT=Fit_OUT(OUT, N.cores=16)

OUT_set=Set_OUT(OUT)

people=Unique_donors(OUT)

phi_matrix=Create_phi_matrix(OUT,people)

save(OUT, file="./OUT/OUT_ALL.NAT.RData")


phi_tmp=reshape2::melt(phi_matrix) %>% filter(value != "NA") %>% select(1,3)

phi_tmp = tibble::column_to_rownames(phi_tmp,var="Var1")

phi=phi_tmp$value; 
names(phi)=gsub("(.*-.*-.*)-.*","\\1",rownames(phi_tmp))

save(phi, file="./OUT/Phi.RData")

ZT=tibble(CaseID=names(phi),phi=round(phi*12/pi,digits = 2))
write.table(ZT,"./OUT/MatchedNAT_sampleTimePoints.txt",sep = "\t",row.names = F,quote = F)




