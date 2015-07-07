
library(kangar00)
library(mboostDevel)

##############################################
##  make GWAS-data object for  simulation   ##
##############################################

library("data.table")  
library("ff")
    
path.sim <- "data/"

      #tex   <- fread(paste(path.sim,"sim2.txt",sep=""), sep='auto', header=FALSE, 
      #         data.table=FALSE)
      #rsnumbers <- read.table("Q:/1_PROJEKTE/sim2.txt",header=FALSE,nrow=1,as.is=TRUE)     
      #rownames in first column, fread doesn't seem to have a "rownames=TRUE" argument
      #rownames(tex) <- tex[,1]   
      #tex           <- tex[,-1]
      #colnames(tex) <- as.character(rsnumbers) #doesn't seem to understand 'header=TRUE'  
      #save(tex, file="tex.rda")
      #select only the part included in simulated pathways:
      #insim <- unique(as.character(anno[,4]))  #129
      #geno <- tex[,colnames(tex)%in%unique(as.character(anno[,4]))] 
      #200 obs. of  129 variables
      #save(geno, file="geno.rda")

anno    <- read.table(paste(path.sim,"sim.6pw.anno.txt",sep=""),
           header=T, as.is=T)
#save(anno, file="anno.rda") 
                     
pheno   <- read.table(paste(path.sim,"sim2.pheno.txt",sep=""),
           header=T, as.is=T)
#save(pheno, file="pheno.rda")
          
load(paste(path.sim,"geno.rda",sep=""))
geno <- as.ffdf(geno)                    
attr(geno,"anno") <- anno

#GWAS object:
sim <- GWASdata(pheno=pheno, geno=geno, desc="sim")





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

### rewire procedure in get_ana
data(hsa04710)
plot(hsa04710)

anno    <- read.table(paste("data/sim.6pw.anno.txt",sep=""),
                      header=T, as.is=T)

table(anno$pathway)

load("data/pathways/hsa00603.RData", verbose = TRUE)
pathway <- p
plot(pathway)

# genes in pathway
net_genes <- get_genes(pathway)
# genes in annotation -> genes with SNPs in GWASdata
anno_genes <- anno$gene[anno$pathway==pathway@id]
# pathway genes that are not in annotation
remov <- which(! net_genes %in% anno_genes)


N <- as.matrix(pathway@adj)
N[N!=0] <- pathway@sign
# include self-interaction
diag(N) <- 1

image(N)

# standard procedure rewire
remov=7
if(length(remov)!=0){
    for(g in remov){
      z <- remov#which(colnames(N)==g)
      vec <- rbind( N[z,],seq(1:length(N[z,])) )
      vec <- data.frame( vec[,vec[1,]!=0])
        #if something must be rewired
        if( length(vec[1,])>1 ){
          for(i in 1:(length(vec[1,])-1) ){
              if((i+1)<=length(vec[1,])){
              for( j in (i+1):length(vec[1,]) ){ #i ist aktuelle edge
               if(N[vec[2,i],vec[2,j]]!=0){print("Edge will be removed!")}
               N[vec[2,i],vec[2,j]] <- N[vec[2,i],vec[2,j]] + vec[1,i]*vec[1,j]
               N[vec[2,j],vec[2,i]] <- N[vec[2,j],vec[2,i]] + vec[1,i]*vec[1,j]}}
         } }
     N <- N[-z,-z]
     N[N>1]    <-  1
     N[N<(-1)] <- -1
     } 
Nold <- N

#' Apply two-step network to rewire network genes if it contains no SNPs in GWASdata (for internal use)
#'
#' @export
#' @author Juliane Manitz
#'
#' @param N adjacency matrix
#' @param remov indication which genes should be removed
#' @return An adjacency matrix containing rewired network
#'
#' @references TODO Newman?
rewire_network <- function(N, remov)
    # early exist if no genes have to be removed
    if(length(remov)=0){ return(N) }

    # identify genes that need to be carried forward to the subnetwork
    ind_sub <- which(N[remov,] != 0)
    # extract the subnetwork
    Nsub <- N[ind_sub, ind_sub]
    # exclude self-interaction
    diag(Nsub) <- 0

    # calculate the two-step network
    Nsub2step <- Nsub %*% Nsub
    Nsub2step[Nsub2step>0] <-  1
    Nsub2step[Nsub2step<0] <- -1

    # check whether interaction types contradict
    check_contradicts <- (Nsub != 0) & (Nsub + Nsub2step == 0)
    if(any(check_contradicts)){
       Nsub2step[(Nsub != 0) & (Nsub + Nsub2step == 0)] <- 0
       message('Interaction types contradict after rewiring: Edges removed.')
    }
    # replace the subnetwork in the adjacency matrix
    N[ind_sub,ind_sub] <- Nsub2step
    # remove the genes and return network
    return(N[-remov,-remov])
}

#' Produce middle part of Network Kernel (for internal use)
#'
#' @export
#' @author Juliane Manitz, Saskia Freytag, Stefanie Freidrichs
#'
#' @param anno \code{data.frame} with annotation information
#' @param SNPset vector with SNPs to be analyzed
#' @param pathway pathway object
#' @return matrix ANA' for inner part of network kernel
get_ana <- function(anno, SNPset, pathway){

    N <- as.matrix(pathway@adj)
    N[N!=0] <- pathway@sign
    if(any(is.na(N))){stop("network information contains missing values")}

    ### remove genes that have no SNPs (not in anno):
    # genes in pathway
    net_genes <- get_genes(pathway)
    # genes in annotation -> genes with SNPs in GWASdata
    anno_sub <- anno[anno$pathway==pathway@id,]
    anno_genes <- unique(anno_sub$gene)
    # pathway genes that are not in annotation
    remov <- which(! net_genes %in% anno_genes)
    # rewire network -> separate function
    N <- rewire_network(N, remov)

    # include selfinteractions for main effects
    diag(N) <- 1
    #ensure positive definiteness by network weighting
    N <- make_posdev(N)

    #A: SNP to gene mapping
    Atab <- table(anno_sub[c('snp','gene')])
    Amat <- as(Atab, 'matrix')
    A <- Amat[SNPset,rownames(N)]    #A is colnames(Z) x rownames(N)

    #A*: size-adjustement for no of SNPs in each gene
    A.star <- A/colSums(A)

    return(A.star %*% N %*% t(A.star))
}






  image(Nold)
  image(Ntmp)


length(get_genes(pathway))
plot(pathway)
Nsub2step
