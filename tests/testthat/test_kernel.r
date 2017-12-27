
library("testthat")

context("Test kernel calculation")

test_that("Test network kernel calculation and dimension", {

  # saved data  
  data(net.kernel.hsa04020)
  # derivation
  data(gwas)
  data(hsa04020)
  net_kernel <- calc_kernel(gwas, hsa04020, knots=NULL, type='net', calculation='cpu')

  expect_equal(dim(net_kernel@kernel)[1],dim(gwas@geno)[1])
  expect_equal(dim(net_kernel@kernel)[2],dim(gwas@geno)[1])

  expect_equal(dim(net_kernel@kernel), net_kernel@kernel)

})

test_that("Test lowrank kernel dimension", {

  data(gwas)
  data(hsa04020)
  square <- calc_kernel(gwas, hsa04020, knots=gwas, type='lin', calculation='cpu')
  dim(square@kernel)
  gwas2 <- new('GWASdata', pheno=pheno[1:10,], geno=geno[1:10,], anno=anno, desc="study 2")
  low_rank <- calc_kernel(gwas, hsa04020, knots = gwas2, type='net', calculation='cpu')
  dim(low_rank@kernel)

  expect_equal(dim(low_rank@kernel)[1],dim(gwas@geno)[1])
  expect_equal(dim(low_rank@kernel)[2],dim(gwas2@geno)[1])
   
})
