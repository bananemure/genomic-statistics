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

###exploring genotype data
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
