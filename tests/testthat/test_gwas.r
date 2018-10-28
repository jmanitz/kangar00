
library("testthat")

context("Test gwas data manipulation")

test_that("Test gwas data manipulation", {

# derivation
data(pheno)
data(geno)
data(anno)
gwas2 <- new('GWASdata', pheno=pheno, geno=geno, anno=anno, desc="some study")
# saved data
data(gwas)

expect_equal(gwas@geno, gwas2@geno)
expect_equal(gwas@anno, gwas2@anno)
expect_equal(gwas@pheno, gwas2@pheno)

})

#test_that("Test SNP and pathway infos", { 
#not stable as pathway info is constantly updated

# derivation
#hsa_info_ex <- pathway_info('hsa04022')
#snp_info_ex <- snp_info(c("rs10243170"))
# saved data
#data(hsa04022_info) 
#data(rs10243170_info) 

#expect_equal(snp_info_ex, rs10243170_info)
#expect_equal(hsa_info_ex, hsa04022_info)

#})


