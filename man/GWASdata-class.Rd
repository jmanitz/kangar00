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
This class creats an object with slots \code{pheno}, \code{geno}, \code{anno} and \desc{desc}. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{GWASdata(...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{pheno}:}{Object of class \code{data.frame} specifying IDs, phenotypes and covariates e.g. ID, pheno, sex, pack.years.}
    \item{\code{geno}:}{Object of class \code{databel} containing Genotype information}
    \item{\code{anno}:}{Object of class \code{data.frame} ontaining annotation file mapping SNP -> genes -> pathways }
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
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GWASdata}}
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("GWASdata")
}
\keyword{classes}
\usage{
GWASdata(...)
}
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
