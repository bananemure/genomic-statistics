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

#show qc criteria computed
names(QC1)

#generate a new data set, which will consist only of people and markers that passed the QC tests
cleanData1<-ge03d2ex[QC1$idok,QC1$snpok]
##nb: The qc can be repeated iteratively until no errors are found
#fix X/Y/mtDNA-errors
cleanData1<-Xfix(cleanData1)
###detach raw phenotypic dataset
detach(phdata(ge03d2ex))
###attach new first clean phenotypic dataset
attach(phdata(cleanData1))

#####################################look at genetic sub-structure ######################
##detect genetic outliers

data1.gkin <- ibs(cleanData1[, autosomal(cleanData1)], weight="freq")
#The numbers below the diagonal show the genomic estimate of kinship 
#(aka 'genomic kinship' or 'genome-wide IBD'), the numbers on the diagonal correspond to 0.5 plus the genomic homozigosity, 
#and the numbers above the diagonal tell how many SNPs were typed successfully for both subjects 
#(thus the IBD estimate is derived using this number of SNPs).

#transform this matrix to a distance matrix using standard R
data1.dist<-as.dist(0.5-data1.gkin)
#perform Classical Multidimensional Scaling:
#By default, the first two principal components are computed and returned.
data1.mds <- cmdscale(data1.dist)
plot(data1.mds)

#identify the points belonging to clusters
km <- kmeans(data1.mds, centers=2, nstart=1000)
cl1 <- names(which(km$cluster==1))
cl2 <- names(which(km$cluster==2))
if (length(cl1) > length(cl2)) {x<-cl2; cl2<-cl1; cl1<-x}

#form a data set which is free from outliers by using only people from the bigger cluster:
cleanData2 <- cleanData1[cl2, ]

#At this stage, we want to allow for HWE checks (we will use only controls and exclude markers with FDR =< 0.2):
QC2 <- check.marker(cleanData2, hweids=(phdata(cleanData2)$dm2==0), fdr=0.2)
summary(QC2)

#Indeed, in the updated data set several markers do not pass our QC criteria and we need to drop a few markers.
cleanData2 <- cleanData2[QC2$idok, QC2$snpok]