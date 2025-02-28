rm(list=ls())
setwd("~/project/CirDatabase_020824/TCGA_NATSample_timePoints_predict/TCGA_NAT_samplePredict/")

source("../Scripts/nconds.R")
source("../Scripts/nconds_functions.R")

library(lmtest)
library(data.table)
library(parallel)
library(combinat)



#Divide E object by Disease status of donors
split_E_Status<-function(E, samp){
  NAT=samp$SAMPID[samp$SAMTYPE=="Solid Tissue Normal"]
  Tumor=samp$SAMPID[samp$SAMTYPE != "Solid Tissue Normal"]
  A=list()
  for(i in names(E)){
    e=E[[i]]
    ee=e$E
    cn=colnames(ee)
    midx=match(NAT, cn)
    midx=midx[!is.na(midx)]
    if (length(midx)>24){
      e$E=ee[,midx]
      A[[paste(i, "NAT", sep="-")]]=e
    }
    fidx=match(Tumor, cn)
    fidx=fidx[!is.na(fidx)]
    if (length(fidx)>24){
      e$E=ee[,fidx]
      A[[paste(i, "Tumor", sep="-")]]=e
    }
  }
  return(A)
}

spliti= function(x, sp, nb){
  v=sapply(strsplit(x,sp),"[[",nb)
  return(v)
}

# Reorganizes the CPM matrix into a list (E), removing unwanted tissues
CPM_to_E<-function(CPM.all.norm, min.samp=20, sep=NULL, samp=NULL){
  E.matrix=list()
  if(!is.null(samp)){
    if(tolower(sep)=="samtype"){
      NAT=samp$SAMPID[samp$SAMTYPE=="Solid Tissue Normal"]
      Tumor=samp$SAMPID[samp$SAMTYPE=="Primary Tumor"]
      for(name in names(CPM.all.norm)){
        E=list()
        E$E=CPM.all.norm[[name]]
        E$E=E$E[,!is.na(E$E[1,])]
        if(!is.null(E$E)){
          E$tissue=name
          ml=fm=list()
          colnames(E$E)=spliti(colnames(E$E), "\\.", 2)
          ml$E=E$E[,colnames(E$E) %in% NAT]
          ml$type="NAT"
          ml$tissue=name
          if(ncol(ml$E)>=min.samp){E.matrix[[paste(name, "NAT", sep="-")]]=ml}
          fm$E=E$E[,colnames(E$E) %in% Tumor]
          fm$type="Tumor"
          fm$tissue=name
          if(ncol(fm$E)>=min.samp){E.matrix[[paste(name, "Tumor", sep="-")]]=fm}
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

#Create the OUT structured file from any CPM and set of phases
Make_big_OUT<-function(E,phi){
  donors=names(phi)
  OUT=list()
  out=list()
  for (i in names(E)) {
    out=E[[i]]
    colnames(out$E)=gsub("\\..*$","", gsub("GTEX.","", colnames(out$E)))
    int_donors=intersect(donors, colnames(out$E))
    if(length(int_donors)>0){
      out$E=out$E[,int_donors]
      out$phi=phi[int_donors]
      OUT[[i]]=out
    }
  }
  return(OUT)
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
#Harmonic regression funcion
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
########

#### Main ####
if(!exists("N.cores")) N.cores = 30 

load("./CPM/NATandTumor.CPM_fullSamplegt25.RData")
CPM.all.norm= get(load("./CPM/NATandTumor.NORM.CPM_fullSamplegt25.RData"))

phi_tmp = read.delim("./OUT/MatchedNAT_sampleTimePointsPirange.txt")
names(phi_tmp)=c("CASEID","Phi")
sample = read.delim("MatchedNAT_Tumor.sampleInfo.txt")

sample_phi = merge(sample,phi_tmp,by="CASEID")

write.table(sample_phi,"./MatchedNAT_Tumor.sampleInfo_withPHIpiranger.txt",sep="\t",quote = F,row.names = F)

phi_tmp=sample_phi %>% select(SAMPID,Phi) %>% tibble::column_to_rownames(.,var="SAMPID")
phi=as.numeric(phi_tmp$Phi)
names(phi)=rownames(phi_tmp)

E=CPM_to_E(CPM.all.norm)
OUT.all=Make_big_OUT(E, phi)
OUT.all=Fit_OUT(OUT.all, N.cores=N.cores)
save(OUT.all, file="./OUT/allSamples.sample.FiT.OUT.RData")

E.status=split_E_Status(E, sample_phi)
OUT.NT=Make_big_OUT(E.status, phi)
OUT.NT=Fit_OUT(OUT.NT, N.cores = N.cores)
save(OUT.NT, file="./OUT/NATvsTumore.sample.FiT.OUT_Status.RData")
load("./OUT/NATvsTumore.sample.FiT.OUT_Status.RData")

qcut=0.2 # BH corrected p-value > 0.2
Rcut=0


tix.c=gsub("-NAT", "",gsub("-Tumor", "", names(OUT.NT)))


SS=list()
for (tx in tix.c[duplicated(tix.c)]){
  MT=subset(sample_phi, PROJID==tx)
  idx=which(tix.c==tx)
  nms=gsub("^.*-", "", names(OUT.NT)[idx])
  nm1=paste(tx,"-NAT", sep="")
  nm2=paste(tx,"-Tumor", sep="")
  
  T1=OUT.NT[[nm1]]$E
  T2=OUT.NT[[nm2]]$E
  P1=OUT.NT[[nm1]]$phi
  P2=OUT.NT[[nm2]]$phi
  S1=match(intersect(MT$SAMPID, colnames(T1)), colnames(T1))
  S2=match(intersect(MT$SAMPID, colnames(T2)), colnames(T2))
  NS=length(S1)
  IRN=intersect(rownames(T1),rownames(T2))
  FM=cbind(T1[IRN,S1], T2[IRN,S2])
  FP=c(P1[S1],P2[S2])
  raw.E=CPM.all.norm[[tx]]
  raw.E=raw.E[IRN,colnames(FM)]
  conds=c(rep("NAT",NS),rep("Tumor",NS))
  ss=nconds(FM,conds=conds,t=FP*12/pi, out.prefix = NULL, N.cores = N.cores)
  
  gn=NULL
  for(nm in c(nm1, nm2)){
    dvt=gsub("^.*-", "", nm)
    out=OUT.NT[[nm]]
    pvals=out$data.fit[,c("qval","amp")]
    rownames(pvals)=out$data.fit[,"genes"]
    pvalus=subset(pvals, qval<qcut & amp>Rcut)
    gn=c(gn, rownames(pvalus))
    coms=intersect(rownames(pvalus), rownames(ss))
    ss$qval=1
    ss$R=0
    ss[coms, c("qval", "R")]=pvals[coms,]
    colnames(ss)[(length(colnames(ss))-1):length(colnames(ss))]=paste(c("qvals", "R"), dvt, sep="_")
  }
  dvt="all"
  out=OUT.all[[tx]]
  pvals=out$data.fit[,c("qval","amp")]
  rownames(pvals)=out$data.fit[,"genes"]
  pvalus=subset(pvals, qval<qcut & amp>Rcut)
  gn=c(gn, rownames(pvalus))
  coms=intersect(rownames(pvals), rownames(ss))
  ss$qval=1
  ss$R=0
  ss[coms, c("qval", "R")]=pvals[coms,]
  colnames(ss)[(length(colnames(ss))-1):length(colnames(ss))]=paste(c("qvals", "R"),dvt, sep="_")
  gnu=unique(gn)
  ss$accepted=0
  ss$accepted[match(gnu,rownames(ss))]=1
  ss$model.c=ss$model*ss$accepted
  ss=ss[,-(ncol(ss)-1)]
  SS[[tx]]=ss
}
save(SS, file="./OUT/SS_SampleTypeModelselect.RData")




table(tem$model)





############################NAT and Tumor seperatedly###

###NAT out fit ###
NAT.fitout_step4 = get(load("./OUT/OUT_ALL.NAT.RData"))
NAT.CPM = get(load("./CPM/NAT.CPM_fullSamplegt25.RData"))
NAT.E=CPM_to_E(NAT.CPM)
NAT.OUT=Make_big_OUT(NAT.E, phi)
NAT.fitout=Fit_OUT(NAT.OUT, N.cores=N.cores)

###Tumore CPM load###
Tumor.CPM = get(load("./CPM/Tumor.CPM_fullSamplegt25.RData"))
Tuomor.E=CPM_to_E(Tumor.CPM)
Tumor.OUT=Make_big_OUT(Tuomor.E, phi)
Tumor.fitout=Fit_OUT(Tumor.OUT, N.cores=N.cores)

####




