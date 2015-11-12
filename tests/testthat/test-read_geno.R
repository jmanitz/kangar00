source("GWASdata.r")
library("testthat")


fread.warn <- "fread caused a warning. Please make sure your file has been read in correctly! \n"
fread.erro <- "fread has stopped due to an error! \n"
fread.load <- "Loading data via fread. If this leads to problems the function will try automatically to load file another method. \n"
read.table.warn <- "read.table caused a warning. Please make sure your file has been read in correctly! \n"
read.table.erro <- "read.table has stopped due to an error! Try to convert your file in a .txt-file and try again. \n"
read.table.load <- "Try reading in file with read.table. Attention: This function is very slow! \n"
read.big.warn <- "read.big.matrix caused a warning. Please make sure your file has been read in correctly! \n"
read.big.erro <- "read.big.matrix has stopped due to an error! \n"
read.big.load <- "Try loading data via read.big.matrix! \n"

context("Test functionality of read_geno:")

test_that("Test for .txt-files without header:", {
  
  expect_is(read_geno, "big.matrix" )
  
  expect_warning
  
  expect_error
  
  
  
})