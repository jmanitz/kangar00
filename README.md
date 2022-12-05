# kangar00

[![Build Status (Linux)](https://travis-ci.org/jmanitz/kangar00.svg?branch=master)](https://travis-ci.org/jmanitz/kangar00)
[![Build status (Windows)](https://ci.appveyor.com/api/projects/status/github/jmanitz/kangar00?branch=master&svg=true)](https://ci.appveyor.com/project/jmanitz/kangar00/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/jmanitz/kangar00/badge.svg?branch=master)](https://coveralls.io/github/jmanitz/kangar00?branch=master)
![GitHub repo
size](https://img.shields.io/github/repo-size/jmanitz/kangar00)
![GitHub issues](https://img.shields.io/github/issues/jmanitz/kangar00)


[![CRAN Status Badge](http://www.r-pkg.org/badges/version/kangar00)](https://CRAN.R-project.org/package=kangar00)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/kangar00)](https://CRAN.R-project.org/package=kangar00)


![CRAN License](https://img.shields.io/cran/l/kangar00) 
![CRAN dependencies status](https://img.shields.io/librariesio/release/CRAN/kangar00)
![Website](https://img.shields.io/website?url=http%3A%2F%2Fkangar00.manitz.org%2F)


Methods to extract information on pathways, genes and SNPs from online databases. It provides functions for data preparation and evaluation of genetic influence on a binary outcome using the logistic kernel machine test (LKMT). Three different kernel functions are offered to analyze genotype information in this variance component test: A linear kernel, a size-adjusted kernel and a network based kernel.

## Citation

To cite the package `kangar00` itself use:

- J. Manitz, S. Friedrichs, P. Burger, B. Hofner, N.T. Ha, S. Freytag, H. Bickeboeller (2017). kangar00: Kernel Approaches for
  Nonlinear Genetic Association Regression, R package version 1.0, https://CRAN.R-project.org/package=kangar00.

The size-adjusted kernel function is introduced in:

- S. Freytag, H. Bickeboeller, C.I. Amos, T. Kneib, M. Schlather (2012). A Novel Kernel for Correcting Size Bias in the Logistic Kernel
  Machine Test with an Application to Rheumatoid Arthritis. Human Heredity, 74, 97-108.

The network-based kernel function is introduced in:

- S. Freytag, J. Manitz, M. Schlather, T. Kneib, C.I. Amos, A. Risch, J. Chang-Claude, J. Heinrich, H. Bickeboeller (2013). A
  Network-Based Kernel Machine Test for the Identifcation of Risk Pathways in Genome-Wide Association Studies. Human Heredity, 76,
  64-75.

The kernel boosting method is introduced in:

- S. Friedrichs, J. Manitz, P. Burger, C.I. Amos, A. Risch, J.C. Chang-Claude, H.E. Wichmann, T. Kneib, H. Bickeboeller, and B. Hofner
  (2017). Pathway-Based Kernel Boosting for the Analysis of Genome-Wide Association Studies. 
  Computational and Mathematical Methods in Medicine, 2017(6742763), 1-17. doi:10.1155/2017/6742763.

Use `toBibtex(citation("kangar00"))` in R to extract BibTeX references.