#source("GWASdata.r")
library("testthat")
library("data.table")

context("Test functionality of read_geno:")

test_that("Test for general errors in read_geno:", {
  
  file.dir.1 <- system.file(package = "kangar00")
  file.dir.2 <- "/inst/data/"
  file.dir   <- paste(file.dir.1, file.dir.2, sep = "")
  
  setwd(file.dir)
  
  file.name <- "no_file.txt"
  file.path <- paste(file.dir, file.name, sep="")
  file.error <- paste("File", file.path, "does not exist!", sep=" ")
  
  expect_error(read_geno(file.path), file.error)
  
  file.name <- "wrong_format.csv"
  file.path <- paste(file.dir, file.name, sep="")
  file.error <- paste("csv as file format is not accepted!")
  
  expect_error(read_geno(file.path), file.error)
  
})

test_that("Test for .txt-files using fread without header and row.names:", {
  
  file.dir.1 <- system.file(package = "kangar00")
  file.dir.2 <- "/inst/data/"
  file.dir   <- paste(file.dir.1, file.dir.2, sep = "")
  
  file.name <- "test_geno_no_header.txt"
  file.path <- paste(file.dir, file.name, sep="")
  
  expect_is(read_geno(file, header = FALSE, row.names = FALSE, use.fread = TRUE), "big.matrix" )
  
  data.vec <- c(0, 1, 0, 2, 1, 0, 1, 2, 0, 1, 1, 0, 2, 2, 1, 2, 2, 1, 0, 1, 1, 0, 2, 0)
  data.ma  <- matrix(data.vec, nrow = 6, ncol = 4)
  data.bma <- as.big.matrix(data.ma)
  
  expec_equal(read_geno(file, header = FALSE, row.names = FALSE, use.fread), data.bma)
  
  msg_warning_ID <- "Your geno file doesn't seem to contain ID numbers. Please make sure that
                   the first row of your data contains ID numbers according to your phenotype file! \n"
  expect_warning(read_geno(file, header = FALSE, row.names = FALSE), msg_warning_ID)
  
})