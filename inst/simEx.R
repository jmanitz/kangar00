library(kangar00)
library(mboostDevel)

path.sim <- "data/"
anno    <- read.table(paste(path.sim,"sim.6pw.anno.txt",sep=""),
                      header=T, as.is=T)
pheno   <- read.table(paste(path.sim,"sim2.pheno.txt",sep=""),
                      header=T, as.is=T)


## check for overlapping snps
tmp <- tapply(anno$snp, anno$pathway, function(x) x)
p <- list()
for (i in 2:6)
    p[[i-1]] <- sapply(tmp[[1]], function(SNP) SNP %in% tmp[[i]])
names(p) <- names(tmp)[-1]
lapply(p, sum)


##### make GWASobjekt for simulated data #####################
if (!file.exists("ff.Rda") || !file.exists("data/sim2.6pw.GWAS.RData")) {
    dat <- databel(paste(path.sim,"sim2",sep=""))
    data <- dat[, unique(anno$snp)]
    data_matrix <- as(data, "matrix")
    ff_dat <- as.ffdf(as.ff(data_matrix, filename = "./small_data.ff"))
    save("ff_dat", file = "ff.Rda")

    load("ff.Rda", verbose = TRUE)
    ## TEST THIS
    # databel2text(dat,file="data_orig.txt")
    ## --> Leerzeichenseparierte Datei
    # data_ffdf <- read.table.ffdf(file = "data_orig.txt", header = TRUE, sep = " ")

    attr(ff_dat, "anno") <- anno

    sim2 <- GWASdata(pheno = pheno, geno = ff_dat,
                     desc = "simulation2.6pathways")
    save(sim2,file="data/sim2.6pw.GWAS.RData")
}

### load prepared data #########################################################
load("ff.Rda", verbose = TRUE)
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

save(hsa00780, file='data/hsa00780.rda')
save()

# GWAS data
sim2
summary(sim2)
GeneSNPsize(sim2)

# pathway information
data(hsa00780)
hsa00780
summary(hsa00780)
plot(hsa00780)

# find a nice/good example
load("data/list.of.pathway.objects.RData", verbose = TRUE)
lapply(p.list[1:100], FUN=analyze)

pdf('ptest.pdf')
lapply(p.list[c(33,64,91,93)],plot) 
dev.off()

hsa04710 <- p.list[[64]]
save(hsa04710, file='../data/hsa04710.rda')

# example
data(hsa04710)
hsa04710
summary(hsa04710)
analyze(hsa04710)
plot(hsa04710)
samp <- sample_genes(hsa04710, no = 3)
plot(hsa04710, highlight.genes = samp)

# sparse matrices for pathways?
require(Matrix)

load('kangar00/data/hsa04710.rda')
object.size(hsa04710)

g <- pathway2igraph(hsa04710)
m <- as_adj(g, type='upper', sparse=TRUE)
object.size(m)

#### now use linear kernel with mboost
data <- sim2@pheno
data$pheno <- as.factor(data$pheno)
data$sex <- as.factor(data$sex)
data$z <- sim2@geno
attr(data$z, "anno")

## not positive definite despite make_posdev
X <- kernel.lin(data$z, pathway = hsa00100)@kernel
e <- eigen(X)
any(e$values < 0)
any(round(e$values, 10) < 0)

## this is positive definite
X <- kernel.lin(data$z, pathway = hsa00603)@kernel
e <- eigen(X)
any(e$values < 0)

mod <- gamboost(pheno ~ bols(sex) + bols(age) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00100, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00603, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00780, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa04122, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00750, args = list()) +
                  bkernel(z, kernel = kernel.lin, pathway = hsa00730, args = list()),
                data = data, family = Binomial())

## selektiert zuerst die beiden informativen pathways!
selected(mod)

extract(mod, "bnames")

## not working
set.seed(1234)
cvr <- cvrisk(mod, grid = 1:200, papply = mclapply)
cvr
plot(cvr)
save("cvr", file = "cvrisk.Rda")
mstop(mod) <- mstop(cvr)

names(coef(mod))


lkmt(pheno ~ sex + age, kernel = kernel.lin(sim2, pathway = hsa00100), GWASdata = sim2)
lkmt(pheno ~ sex + age, kernel = kernel.lin(sim2, pathway = hsa00603), GWASdata = sim2)
lkmt(pheno ~ sex + age, kernel = kernel.lin(sim2, pathway = hsa00780), GWASdata = sim2)
lkmt(pheno ~ sex + age, kernel = kernel.lin(sim2, pathway = hsa04122), GWASdata = sim2)
lkmt(pheno ~ sex + age, kernel = kernel.lin(sim2, pathway = hsa00750), GWASdata = sim2)
lkmt(pheno ~ sex + age, kernel = kernel.lin(sim2, pathway = hsa00730), GWASdata = sim2)


### test gputools


require(gputools)

# test matrix multiplication
matA <- matrix(runif(5000*3000),5000,3000)
matB <- matrix(runif(5000*3000),3000,5000)

system.time(cpuMatMult(matA,matB))
system.time(gpuMatMult(matA,matB))

# linear kernel
zi <- sample(c(0,1,2),size=4000*500000))
Z <- matrix(zi, 4000,500000)
system.time(gpuMatMult(Z,t(Z)))
system.time(cpuMatMult(Z,t(Z)))

