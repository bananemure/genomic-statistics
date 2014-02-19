rm(list=ls())
setwd('genomics_biostats/GWAS/')
require(BiocInstaller)
#biocLite('BiocUpdate')
#biocLite('gage')
#biocLite('GenABEL')


########################################
#                                     ##
#             gwas-GenABEL            ##
#                                     ##
########################################
library(GenABEL)
#load data
data(srdta)
#look at the structure of srdata genabel type data
str(srdta)
####explore phenotype data
head(phdata(srdta))
names(phdata(srdta))

###explore genotype data
str(gtdata(srdta))
#number of participants
nids(gtdata(srdta))
#number of SNPs
nsnps(gtdata(srdta))
#ids of individuals
idnames(srdta)
#sex of individuals
male(srdta)
#SNPs names
snpnames(srdta)
#chromosome of bearing SNPs
chromosome(srdta)
#map position of SNPs
map(srdta)
#coding of SNPs
coding(srdta)
#the strand on which the coding is reported
strand(gtdata(srdta)): # '-','+',or 'u' for missing
#ref allele of SNPs
refallele(srdta)
#effective allele of SNPs
effallele(srdta)

#####acessing and modifying phenotype data
#exemple adding a new colum of square of age
age2<-srdta@phdata$age^2
srdta<-add.phdata(srdta,newphdata=age2,name='squaredAge')
names(phdata(srdta))
# delete the squareAge variable
srdta<-del.phdata(srdta,what='squaredAge',all=F)
names(phdata(srdta))

###################exploring genotype data##################################

#view genabel data as genotype like data: same as genetics package(obselete) 
head(as.character(srdta),2)
#translate GenABEL-package genetic data to the format used by "genetics" library
geneticsLikeFormat<-as.genotype(gtdata(srdta))
#translate GenABEL-package data to the format used by "haplo.stats"
haploLikeFormat<-as.hsgeno(gtdata(srdta))

#####MINING srdta genotype data
#see actual genotype of people p141,p147,p2000 and their SNPs rs10,rs29
as.character(srdta[c("p141", "p147", "p2000"), c("rs10", "rs29")])
### exemple: Exploring rs114:
#What is the coding and which allele is the reference one?
c(coding(srdta[,'rs114']),refallele(srdta[,'rs114']))
#What is the frequency of non-reference ("effect") allele in total sample?
#check first the mapping of character coding with numeric coding: homozygote effalle are mapped to 2 and heterozygote are mapped to 1
#homozygote refallele are mapped to 0 ==> freq=homozygote+1/2(heterozygote)==> sum/2*population because homozygote effalle are counted in double
table(as.character(gtdata(srdta[,'rs114'])),as.numeric(gtdata(srdta[,'rs114'])))

effQty<-as.numeric(gtdata(srdta[,'rs114']))
freqEffect<-mean(effQty,na.rm=T)/2
#compute the frequency of the effect allele of "rs114" in males
effQtyMale<-as.numeric(gtdata(srdta[male(srdta)==1,'rs114']))
freqEffectMale<-mean(effQtyMale,na.rm=T)/2
#compute the same statistic for reference allele
refQty<-1-effQty
refQtyMale<-1-effQtyMale

#### To test for HWE in first 10 SNPs in total sample: check the Pexact column
summary( gtdata(srdta[, 1:10]) )

#######exemple summaries for individuals older than 60
defAge<-which(srdta@phdata$age>60)
### To test for HWE in first 10 SNPs in defAge people: check the Pexact column
summary( gtdata(srdta[defAge, 1:10]) )

##### Look for call rate in whole sample
summaryWhole<-summary(gtdata(srdta))
names(summaryWhole)
calRateWhole<-summaryWhole$CallRate
hist(calRateWhole)
# how many markers had call rate lower than, say, 93%, between 93 and 95%, between 95 and 99% and more than 99%
catable(calRateWhole)
### how many deviate from HWE
devHWE<-summaryWhole$Pexact
beferonniCorrHWE<-0.05/nsnps(srdta)  #beferroni correction for multi-testing
#table to count with respect to Pvalue
catable(devHWE,c(beferonniCorrHWE,.01,.05,.1))
#cumulate table
catable(devHWE,c(beferonniCorrHWE,.01,.05,.1),cumulative=T)

###investigate the minor allele frequency distribution
alFreq<-summaryWhole$Q.2
malFreq<-pmin(alFreq,1. - alFreq)
hist(malFreq)
catable(alFreq, c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99))
catable(malFreq, c(0, 0.01, 0.05, 0.1, 0.2), cum=TRUE)
#visualize the distribution
par(mfrow=c(1,2))
hist(alFreq)
hist(malFreq)

##another type of summary: per id summary
#let us look at the distribution of heterozygosity for whole sample
het<-perid.summary(gtdata(srdta))$Het
summary(het)
catable(het, c(0.1, 0.25, 0.3, 0.35, 0.5))
par(mfrow=c(1,1))
hist(het)




