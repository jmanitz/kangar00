# Step 0: Load different parts of code

source('Q:/2015/kangar00/R/pathway.r')
source('Q:/2015/kangar00/R/kernel.r')
source('Q:/2015/kangar00/R/lkmt.r')
source('Q:/2015/kangar00/R/GWASdata.r')

# If package is build use:
# load(paste(system.file(package = "kangar00"), "/R/pathway.r", sep=""))
# load(paste(system.file(package = "kangar00"), "/R/kernel.r", sep=""))
# load(paste(system.file(package = "kangar00"), "/R/lkmt.r", sep=""))
# load(paste(system.file(package = "kangar00"), "/R/GWASdata.r", sep=""))


# Step 1: Load test data 

load("Q:/2015/kangar00/data/geno.rda")

load("Q:/2015/kangar00/data/anno.rda")
  
load("Q:/2015/kangar00/data/pheno.rda")

# If package is build use:
# load(paste(system.file(package = "kangar00"), "/data/geno.rda", sep=""))
# load(paste(system.file(package = "kangar00"), "/data/anno.rda", sep=""))
# load(paste(system.file(package = "kangar00"), "/data/pheno.rda", sep=""))

# Step 2: Create sample GWAS object

desc <- "This is an fictive sample for a GWAS object"  

gwas2 <- new("GWASdata")
gwas2@geno <- geno
gwas2@anno <- anno
gwas2@pheno <- pheno
gwas2@desc <- desc

gwas <- GWASdata(geno = geno, anno = anno, pheno = pheno, desc = desc)

# Method 3 does not work for reasons
#gwas3 <- GWASdata(geno, anno, pheno, desc)

# Step 3: Calculate p-vaules
pws <- unique(anno$pathway)

plin <- numeric()
psia <- numeric()
pnet <- numeric()

for(i in 1:length(pws)){
  p <- get(load(paste("Q:/2015/Island/pathways/", pws[i], ".RData", sep = "")))
  # p <- get(load(paste(system.file(package = "kangar00"), "/inst/data/pathways/", pws[i], ".RData", sep = "")))
 
  KL <- lin_kernel(gwas, p, parallel = 'none')
  lobj <- lkmt(pheno ~ 1+age+sex, KL, gwas)
  plin[i] <- lobj@p.value
  
  KS <- sia_kernel(gwas,  p, parallel = 'none')
  sobj <- lkmt(pheno ~ 1+age+sex, KS, gwas)
  psia[i] <- sobj@p.value
  
  KN <- net_kernel(gwas, p, parallel = 'none')
  nobj <- lkmt(pheno ~ 1+age+sex, KN, gwas)
  pnet[i] <- nobj@p.value
  
}

