#############################################################
#############################################################                          
######            GWAS ANALYSIS II WITH GenABEL       #######
#############################################################
#############################################################
#This analysis will use a stratified population dataset All these populations originate from the same base population some
#generations ago. Some of these populations mantained large size and some were small. 
#There was little (2.5%) migration between populations. #wo traits ('quat' and 'bint') ara available for analysis. 
#Our goal is to Investigate relations between phenotypes and covariates. We will Perform association analysis.
rm(list=ls())
library(GenABEL)
load('RData/stratified.RData')

#small exploration of the data
names(phdata(strdat))

#general descriptive stats for phenotypic data
descriptives.trait(strdat)

##################################################################################
#What covariates are significantly associated with the traits?
summary(lm(quat~sex+age+age2,data=phdata(strdat)))      #quantitative trait
summary(glm(bint~sex+age+age2,data=phdata(strdat),family=binomial))   #qualitative trait

#How many SNPs and IDs are presented in the data set
nsnps(strdat)
nids(strdat)

#How many SNPs and IDs pass the quality control at snp and id call rate of 0.98
qcStrat<-check.marker(strdat,callrate=.98,perid.call=.98,p.level=0)
strData1<-strdat[qcStrat$idok,qcStrat$snpok]
attach(strData1@phdata)
nsnps(strData1)
nids(strData1)

#check for evidence of stratification
descriptives.marker(strData1)[2]

#evidence that the test statistics for trait quat is inflated
qts.quat<-qtscore(quat~sex+age+age2,data=strData1)
lambda(qts.quat)

#evidence that the test statistics for trait bint is inflated
qts.bint<-qtscore(bint~sex+age+age2,data=strData1)
lambda(qts.bint)

#How many genetically distinct populations are present in the dataset
gkin.str<-ibs(data=strData1,weight='freq')
dist.str<-as.dist(.5-gkin.str)
pcs.str<-cmdscale(dist.str,k=10)
plot(pcs.str)

