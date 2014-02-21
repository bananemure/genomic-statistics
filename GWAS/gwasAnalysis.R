#############################################################
#############################################################                          
######            GWAS ANALYSIS WITH GenABEL          #######
#############################################################
#############################################################
rm(list=ls())
library(GenABEL)
data(ge03d2ex)
#small exploration of the data
names(phdata(ge03d2ex))

##############descriptive statisitcs##########################
attach(phdata(ge03d2ex))
#general descriptive stats for phenotypic data
descriptives.trait(ge03d2ex)
#by diabete type2
descriptives.trait(ge03d2ex,by=dm2)
##descriptive stats for the marker (genotypic data)
descriptives.marker(ge03d2ex)
#filter by control: dm==0
descriptives.marker(ge03d2ex,idsubset=(dm2==0))
#let us look at HWE in cases: dm2==1
descriptives.marker(ge03d2ex,idsubset=(dm2==1))[2]
#let us look at HWE in controls: dm2==0
descriptives.marker(ge03d2ex,idsubset=(dm2==0))[2]

##GWA scan using the raw (before quality control) data: fast method
rawGwaStats<-qtscore(dm2,ge03d2ex,trait.type='binomial') #univariate analyis: dm2

#checking for inflation without plot
lambda(rawGwaStats)
# checking for inflation with plot  
estlambda(rawGwaStats@results$P1df,plot=T)

##manhattan plot of the results
plot(rawGwaStats)

##add corrected P-values to the plot
add.plot(rawGwaStats, df="Pc1df", col=c("lightblue", "lightgreen"))

##show result of the best 10 P-value(by default=10)
descriptives.scan(rawGwaStats,sortby='P1df')

########### Quality Control ##################
#because of possible inflation or stratification or data error, QC Analysis is important first
## first qc analysis with very low P-value=0: do not check for HWE
#this function allow to set a lot of arguments specifically the platform you use of genotyping: this example is very basic
QC1<-check.marker(ge03d2ex,p.level=0)
#summary of the QC by marker and by person: fails statistics
summary(QC1)
<<<<<<< HEAD
#show qc criteria computed
names(QC1)

#generate a new data set, which will consist only of people and markers that passed the QC tests
cleanData1<-ge03d2ex[QC1$idok,QC1$snpok]
