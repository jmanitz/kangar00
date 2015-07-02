\name{GWASdata-class}
\Rdversion{1.1}
\docType{class}
\alias{GWASdata-class}
\alias{GWASdata}
\alias{GeneSNPsize,GWASdata-method}
\alias{show,GWASdata-method}
\alias{summary,GWASdata-method}

\title{Class \code{GWASdata}}
\description{
This class creats an object with slots \code{pheno}, \code{geno} and \code{desc}. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{GWASdata(pheno, geno, desc)}.
}
\section{Slots}{
  \describe{
    \item{\code{pheno}:}{Object of class \code{data.frame} specifying IDs, phenotypes and covariates e.g. ID, pheno, sex, pack.years.}
    \item{\code{geno}:}{Object of class \code{databel} containing Genotype information and annotation file.}
    \item{\code{desc}:}{Object of class \code{character} giving GWAS description and information }
  }
}
\section{Methods}{
  \describe{
    \item{GeneSNPsize}{\code{signature(object = "GWASdata")}: Create list of pathway names and SNP size and genes size i.e. counts number of Genes/SPNs in each pathway}
    \item{show}{\code{signature(object = "GWASdata")}: Summarize GWASdata information. }
    \item{summary}{\code{signature(object = "GWASdata")}: ... }
	 }
}

\seealso{
\code{\link{GWASdata}}
}
\examples{
showClass("GWASdata")
}
\keyword{classes}
\usage{
GWASdata(pheno, geno, desc) 
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
}
