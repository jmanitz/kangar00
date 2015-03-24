library(kangar00)
library(mboostDevel)

##### make GWASobjekt for simulated data #####################
path.sim <- "data/"
anno    <- read.table(paste(path.sim,"sim.6pw.anno.txt",sep=""), header=T, as.is=T)
pheno   <- read.table(paste(path.sim,"sim2.pheno.txt",sep=""), header=T, as.is=T)
dat     <- databel(paste(path.sim,"sim2",sep=""))
attr(dat, "anno") <- anno

sim2 <- GWASdata(pheno = pheno, geno = dat,
                 desc = "simulation2.6pathways")
save(sim2,file="data/sim2.6pw.GWAS.RData")


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

#### now use linear kernel with mboost

## SNPset <-sim@anno$rsNumber[which(sim@anno$Pathway==pnet@id)]
## z <- as(as(sim@geno[,as.character(SNPset)],'matrix'), "Matrix")

data <- sim2@pheno

data$pheno <- as.factor(data$pheno)
data$sex <- as.factor(data$sex)
data$z <- sim2@geno
attr(data$z, "anno")

mod <- gamboost(pheno ~ sex + age +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00100, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00603, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa04122, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00750, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00730, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00780, args = list()),
                data = data, family = Binomial())

