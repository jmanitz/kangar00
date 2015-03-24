
library(DatABEL) # genotype data objects
library(Matrix)
library(lattice) # levelplot
library(kangar00)
library(mboostDevel)

##### make GWASobjekt for simulated data #####################
path.sim <- "..."
anno    <- read.table(paste(path.sim,"sim.6pw.anno.txt",sep=""), header=T, as.is=T)
pheno   <- read.table(paste(path.sim,"sim2.pheno.txt",sep=""), header=T, as.is=T)
dat     <- databel(paste(path.sim,"sim2",sep=""))
attr(dat, "anno") <- anno

sim2 <- GWASdata(pheno=pheno, geno=dat, desc="simulation2.6pathways")
save(sim2,file="sim2.6pw.GWAS.RData")


### load prepared data #########################################################
load("data/sim2.6pw.GWAS.RData", verbose = TRUE)
load("data/pathways/hsa00100.RData", verbose = TRUE)
hsa00100 <- p
load("data/pathways/hsa00603.RData", verbose = TRUE)
hsa00603 <- p
load("data/pathways/hsa04122.RData", verbose = TRUE)
hsa04122 <- p
load("data/pathways/hsa00750.RData", verbose = TRUE)
hsa00750 <- p
load("data/pathways/hsa00730.RData", verbose = TRUE)
hsa00730 <- p
load("data/pathways/hsa00780.RData", verbose = TRUE)
hsa00780 <- p

# GWAS data
sim2
summary(sim2)
GeneSNPsize(sim2)

# pathway information
hsa00780
summary(hsa00780)
plot(hsa00780)

## # linear kernel
## kern <- kernel.lin(sim, pnet)
## plot(kern)
## kern
## summary(kern)
##
## # size-adjusted
## kern <- kernel.sia(sim,pnet)
## kern
## summary(kern)
## plot(kern)
##
## # network kernel
## kern <- kernel.net(sim,pnet)
## summary(kern)
## plot(kern)
## plot(kern,hclust=TRUE)
##
## # define LKMTtest
## source('kangar00/R/lkmt.r')
## m <- lkmt(pheno ~ 1+sex+age, kern, sim)
##
## m
## summary(m)
##
##
## #### now use linear kernel with mboost
## library("mboostDevel")
##
## z <- sim@geno
## attr(z, "anno") <- sim@anno
##
## ## SNPset <-sim@anno$rsNumber[which(sim@anno$Pathway==pnet@id)]
## ## z <- as(as(sim@geno[,as.character(SNPset)],'matrix'), "Matrix")
##
## data$pheno <- as.factor(data$pheno)
## data$sex <- as.factor(data$sex)
## data$z <- z
## attr(data$z, "anno")
##
## # ND <- data[1:5, ]
## # str(ND)
##
## mod <- gamboost(pheno ~ bkern(z, kernel = kernel.lin, pathway = pnet, args =
##                                   list()) + sex + age,
##                 data = data, family = Binomial())
## plot(mod)
