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

###detach the old and attach the new phenotype data in the search path
detach(phdata(cleanData1))
attach(phdata(cleanData2))

#Before proceeding to GWA, let us check if complete QC improved the fit of genetic data to HWE:
#you can see that the problems with HWE are apparently fixed; we may guess
#that these were caused by Wahlund's effect.
descriptives.marker(cleanData2)[2]

descriptives.marker(cleanData2[dm2==1, ])[2]
descriptives.marker(cleanData2[dm2==0, ])[2]


############################ --- GWAS ANALYSIS --- ###########################
################################**************################################
#Let us start again with descriptives of the phenotypic and marker data: look at Ptt(t.test)
descriptives.trait(cleanData2,by.var=dm2)

###score test on cleaned data2
data2.qt <- qtscore(dm2, cleanData2, trait="binomial")
#check lambda: there is still some inflation, which is explained by the fact that we investigate
#only a few short chromosomes with high LD and few causative variants
lambda(data2.qt)
#Produce the association analysis plot
plot(data2.qt, df="Pc1df")

##scan summary to show the top 10:
descriptives.scan(data2.qt,top=10,sortby='Pc1df')

##compute empirical GW significance for diabete(dm2)
data2.qte<-qtscore(dm2,trait.type='binomial',times=1000,quiet=T,data=cleanData2)
descriptives.scan(data2.qte,sortby='Pc1df')

##adjust the GW significance for dm2 ~ sex and age:
data2.qtsa<-qtscore(dm2~sex+age,data=cleanData2,trait.type='binomial',times=1000,quiet=T)
descriptives.scan(data2.qtsa,top=10,sortby='Pc1df')

##Finally, let us do stratified (by BMI) analysis. We will contracts obese (BMI >30) cases to all controls.
data2.qtsaBmi<-qtscore(dm2~sex+age,data=cleanData2,idsubset=((bmi>=30 & dm2==1) |dm2==0),trait.type='binomial',times=1000,quiet=T)
descriptives.scan(data2.qtsaBmi,top=10,sortby='Pc1df')



#########################################################
####        Genetics stratification             #########
#########################################################
##############WE USE CleanData1 for this analysis #######
#########################################################

####Part1: ETHNIC ADMIXTURE########################################################################################
#detach cleanData2 and attach cleanData1
detach(phdata(cleanData2))
attach(phdata(cleanData1))

# fist check for inflation
data1.qt<-qtscore(dm2,cleanData1,trait.type='binomial')
lambda(data1.qt)

#####1- structured analysis
#define which populations (0 or 1) are in cluster1 (cl1)
pop<-as.numeric(idnames(cleanData1) %in% cl1)
#look how pop is distributed among cases/control
table(pop,dm2,useNA='ifany')
#structure association test
data1.sa <- qtscore(dm2,data=cleanData1,strata=pop,trait.type='binomial')
lambda(data1.sa)
### compare this result with results where ouliers(cl1) were removed 
par(mfcol=c(1,3))
plot(data1.qt,ylim=c(1,6))
plot(data2.qt,ylim=c(1,6))
plot(data1.sa,ylim=c(1,6))
par(mfcol=c(1,1))
#In this case, there is little difference, because all people belonging to the smaller sub-population are cases.

#Other way to adjust for genetic (sub)structure is to apply the method of Price et al. (EIGENSTRAT), which make 
#use of principal components of the genomic kinship matrix to adjust both phenotypes and genotypes for possible stratification

data1.eg <- egscore(dm2,data=cleanData1,kin=data1.gkin)
lambda(data1.eg)
plot(data1.eg,ylim=c(1,6))

#Now let us apply adjustment for the stratification by use of the principal
#components of genetic variation. For that we first need to extract the principal
#components of genetic variation by constructing the distance matrix
dst <- as.dist(0.5-data1.gkin)
#performing the classical multidimensional scaling
pcs <- cmdscale(dst,k=10)
#Now we can use these PCs for adjustment:
data1.pca <- qtscore(dm2~pcs[,1]+pcs[,2]+pcs[,3],cleanData1,quiet=T)
lambda(data1.pca)
plot(data1.pca)

##Finally, let us use the full genomic kinship matrix for the adjustemnt for
##population structure. First, let us estimate the polygenic model with
h2a <- polygenic(dm2,kinship.matrix=data1.gkin,data=cleanData1)
#The resulting 'heritability' estimate is
h2a$esth2

####Now we can perform MIXED polygenic MODEL (Y = mu + G + e) approximation analysis using mmscore function
data1.mm <- mmscore(h2a,cleanData1)
lambda(data1.mm)
###plot all three methods together
par(mfcol=c(1,3))
plot(data1.eg,ylim=c(1,6))
plot(data1.pca,ylim=c(1,6))
plot(data1.mm,ylim=c(1,6))
par(mfcol=c(1,1))

#############   PART 2: Analysis of family data #####################################################################
#we use a new data set: first we clean everything
rm(list=ls())
data(ge03d2.clean)
#subset data
erfs <- ge03d2.clean[1:100,]
#add a new attribute to erfs data
erfs<-add.phdata(erfs,phdata(erfs)$weight,'qtbas') #or phdata(erfs)$qtbas <- phdata(erfs)$weight
#simulate matrix of pedegree
pkins <- matrix(rnorm(nids(erfs)^2,sd=0.01),ncol=nids(erfs),nrow=nids(erfs))

#basic data exploratory
nids(erfs)
class(pkins)
class(erfs)
nsnps(erfs)
#check the distribution of SNPs by chromosome:
table(chromosome(erfs))
#descriptive stats
descriptives.marker(gtdata(erfs))
#histogram of the distribution of the pedegree kinship coefficients
hist(pkins[lower.tri(pkins)])

###Now estimate genomic kinship matrix using autosomal data
gkins<-ibs(erfs[,autosomal(erfs)],weight='freq')
#summary of genomic kinship
summary(gkins[lower.tri(gkins)])
hist(gkins[lower.tri(gkins)],breaks=100)

#relations between genomic and pedigree kinship
plot(pkins[lower.tri(pkins)],gkins[lower.tri(gkins)])
cor(pkins[lower.tri(pkins)],gkins[lower.tri(gkins)])

####analyse the data using plain GC method:
qts <- qtscore(qtbas,data=erfs)
lambda(qts)$est
#The top 10 hits from GWA analysis
descriptives.scan(qts,sortby='Pc1df')
#estimate genome-wide empirical significance with times argument which tells the number of permutations
qts.e <- qtscore(qtbas,data=erfs,times=200,quiet=TRUE)
descriptives.scan(qts.e,top=15,sortby='Pc1df')

##Let us estimate polygenic model
h2 <- polygenic(qtbas,kin=gkins,data=erfs)
#The results of estimation are contained in "h2an" element of the resulting object
#In the "estimate" list, the MLEs shown correspond to intercept, heritability and total variance
h2$h2an

#Let us run FASTA test using estimated polygenic model, as specified by h2
mms <- mmscore(h2,data=erfs)
lambda(mms)$est 
descriptives.scan(mms,sort="Pc1df")

#However, we can not estimate genome-wide significance with FASTA, because the data structure is not exchangeable.
#Using GRAMMAS method, you can estimate nominal P-values
grs <- qtscore(h2$pgres,data=erfs,clam=FALSE)
lambda(grs)$est
