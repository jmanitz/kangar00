\name{GWASdata}
\alias{GWASdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
New GWASdata object
}
\description{
Creates new object of class \code{GWASdata} with slots \code{pheno}, \code{geno} and \code{desc}. 
}
\usage{
GWASdata(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
  Data.frame specifying IDs, phenotypes and covariates e.g. ID, pheno, sex, pack.years.
}
  \item{geno}{
  Genotype information in databel format.
}
  \item{anno}{
  Data.frame containing annotation file mapping SNP -> genes -> pathways
  }
  \item{desc}{
  Character giving GWAS description and information
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}

\keyword{GWASdata}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
