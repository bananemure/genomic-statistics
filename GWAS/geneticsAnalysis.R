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

unlink("RData",recursive=TRUE,force=TRUE)
dir.create("RData")

myDownloads <- function(baseUrl,baseLocal,files) {
  for (cFile in files) {
    cFileUrl <- paste(baseUrl,cFile,sep="")
    cFileLocal <- paste(baseLocal,cFile,sep="")
    tryDownload <- try(
      download.file(url=cFileUrl,destfile=cFileLocal)
    )
    if ( is(tryDownload,"try-error") )
      stop(paste("can not download",cFileUrl,"into",cFileLocal,":",tryDownload))
  }
}

baseUrl <- "http://www.genabel.org/sites/default/files/data/"
baseLocal <- "RData/"
dataFiles <- c(
  "assocbase.RData",
  "popdat.RData",
  "mach1.out.mlinfo",
  "mach1.mldose.fvi",
  "mach1.mldose.fvd",
  "rcT.PHE",
  "gen0.illu",
  "gen0.illuwos",
  "gen0.tped",
  "gen0.tfam",
  "gen0.ped",
  "map0.dat",
  "emap0.dat",
  "phe0.dat",
  "ImputedDataAnalysis.RData")
myDownloads(baseUrl,baseLocal,dataFiles)


#######################################################
#   simple genetics analysis with package genetics  ###
#######################################################
#genetics package is outdate: a new alternative is 'HardyWeinberg' package

library('genetics') 
help(package='genetics')
load("RData//assocbase.RData")
View(assoc)

#convert snps into proper 'genetics' format

assocg <- assoc
assocg$snp4 <- as.genotype(assocg$snp4)
assocg$snp5 <- as.genotype(assocg$snp5)
assocg$snp6 <- as.genotype(assocg$snp6)

#explore assoc data 
class(assocg)
names(assocg)
dim(assocg)
class(assoc$snp4)
attach(assocg)
summary(snp4)
#Let us check if HWE holds for the SNPs described in this data frame
HWE.exact(snp4)
#Let us check if the there is LD between snp4 and snp5:
LD(snp4,snp5)

#######Genetic association analysis
#relation between quantitive trait(qt) and SNPs
mg <- lm(qt~snp4)
summary(mg)

#model with sex covariate
summary(lm(qt~snp4+sex))
#we add interaction between sex and snp4
summary(lm(qt~snp4*sex))

###Test different models: additive,dominant and recessive
##Model additive: the mean value of the trait for heterozygous genotypes is right in between the two homozygotes
#first we convert into numeric [recoding]
add4 <- as.numeric(snp4)-1
#check the quality of recoding
table(snp4,add4,useNA='ifany') 
#now we test the additive model 
summary(lm(qt~add4))

##Model dominant:the means of genotypes 'AA' and 'AB' are the same. recessif for BB
#recoding
dom4<-as.numeric(snp4)
dom4<-replace(dom4,which(snp4=='A/A'|snp4=='A/B'),1)
dom4<-replace(dom4,which(snp4=='B/B'),0)
#check the quality of recoding
table(snp4,dom4,useNA='ifany') 
#test the dominant model
summary(lm(qt~dom4))


####Association with a binary outcome ex:aff(1/0)
summary(glm(aff~snp4,family="binomial"))
#a test of global significance of the SNP effect
anova(glm(aff~snp4,family="binomial"),test="Chisq")

#############Exploration part 2 ##########################################
rm(list=ls())
load('RData//popdat.RData')
attach(popdat)
#compute number of snp in popdat
length(grep('snp',names(popdat)))
#the frequency (proportion) of snp1 allele 'A'
summary(popdat$snp1)
#frequency of 'A' in affected (aff==1)
affected<-which(popdat$aff==1)
summary(popdat$snp1[affected])
#How many cases and controls are present in the data set
table(popdat$aff)

#test for significant association between affection and sex
t<-table(aff,sex)
fisher.test(t)
summary(glm(aff~sex,family=binomial()))

#is association between the disease and qt significant?
summary(glm(aff~qt,family=binomial()))

# test significant association beteween SNPs and affection
for (i in 1:10) {
  snpname <- paste("snp",i,sep="")
  cat ("\n       *******",toupper(snpname),"******      ")
  cat("\nTesting association between aff and SNP",snpname,":\n -------------------------------------- \n")
  print(anova(glm(aff~get(snpname),family=binomial),test="Chisq"))
  cat ("\n============================================\n\n")
  #print(summary(lm(qt~get(snpname)))$coef)
}


