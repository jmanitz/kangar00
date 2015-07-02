\name{GWASdata}
\alias{GWASdata}
\title{
New GWASdata object
}
\description{
Creates new object of class \code{GWASdata} with slots \code{pheno}, \code{geno} and \code{desc}. 
}
\usage{
GWASdata(pheno, geno, desc) {
}
\arguments{
  \item{snps}{
  data.frame specifying ids, phenotypes and covariates e.g. ID, pheno, sex, pack.years.
}
 \item{geno}{
  genotype information in ffdf format. Needs an attribute data frame 'anno' including the columns 'pathway' (pathway id in hsa number as used in the KEGG database), 'gene' (gene name), 'chr' (number of chromosome), 'snp' (rs identifier of genetic marker), 'position' (the markers position (basepairs)).     
}
  \item{desc}{
  character giving GWAS description and information.
}
}

\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner
}

\seealso{
For class information see \code{\link{GWASdata-class}}.
}
\keyword{GWASdata}

