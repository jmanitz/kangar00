
library("testthat")

context("Test LKMT calculation")

test_that("Test LKMT calculation: statistic, df and p-values", {

  # saved data  
  data(lkmt.net.kernel.hsa04020)
  # derivation
  data(hsa04020)
  data(gwas)
  # compute kernel
  net_kernel <- calc_kernel(gwas, hsa04020, knots=NULL, type='net', calculation='cpu')
  # perform LKMT 
  res <- lkmt_test(pheno ~ sex + age, net_kernel, gwas, method='satt')

  expect_equal(lkmt.net.kernel.hsa04020@statistic, res@statistic)
  expect_equal(lkmt.net.kernel.hsa04020@df, res@df)
  expect_equal(lkmt.net.kernel.hsa04020@p.value, res@p.value)

})


