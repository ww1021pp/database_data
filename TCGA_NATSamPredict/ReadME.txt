In the directory, it contained the preocess to predict the TCGA NAT samples time points to do rhythm comparative 
between Cancer and matched NAT sample to get the genes loss rhythm in disease status.

### Step1 ####
step1: Download TCGA raw counts for all Project ID.
DownloadRawCountsFromTCGA.R

###Step2 ####
step2: As we need do time predict, if the NAT samples number is too small, then may be not work well, so we filter Project ID with at least 25 NAT samples.
TCGA_sampleStastics.R 

###Step3 ####
Step3: use EdgeR to get the CPM for all matched NAT and tumor samples, and Remove covariates using a linear regression with age, sex
TCGANATSamplegt25rawCountsTOTPM.R

###Step4 ####
Step4: Use CHIRAL(https://github.com/naef-lab/CHIRAL/tree/master/Pkg/CHIRAL) a Bayesian method to infer the circular coordinates of a set samples.
CHIRAL_NATsamplepredict.R

###Step5 ####
Step5: After process Step4, we can get the circular coordinate of TCGA NAT samples(number of samples gt 25), 
we can use DryR to an R package that provides the statistical framework to assess differential rhythmicity of 
a time series of RNA-Seq data with two and more conditions

Infinally, we just use the model select to compare the rhythmic different between NAT and Cancer Condition.

####Step6 #### plot the gene expression based on Modelselect results##
Step6.1 based on the facet_wrap 

 



