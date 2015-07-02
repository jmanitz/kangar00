\name{kernel.lin}
\alias{kernel.lin}
\title{
Linear Kernel
}
\description{
Creates the linear kernel to be evaluated in the logistic kernel machine test. 
}
\usage{
kernel.lin(data, pathway, parallel = c("none", "cpu", "gpu"))
}
\arguments{
  \item{data}{
  Object of the class \code{GWASdata} containing the genotypes of the individuals for which a kernel will be calculated.  
}
  \item{pathway}{
  Object of the class \code{pathway}containing the genotypes of the individuals for which a kernel will be calculated.
}
  \item{parallel}{
  Argument specifying whether to calculate kernels on several CPUs in parallel (\code{"cpu"}) or on ond CPU (\code{"none"}) or on a GPU (\code{"gpu"}). 
}
}
\value{
Returns object of class kernel, including the similaritymatrix of the corresponding pathway for the considered individuals.  
}
\references{
Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin X: Powerful SNP-Set Analysis for Case-Control Genome-Wide Association Studies. Am J Hum Genet 2010, 86:929-42
}
\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner with contributions by Saskia Freytag and Ngoc-Thuy Ha.
}

\seealso{
\code{\link{GWASdata}}, \code{\link{pathway}}
}

\keyword{LinearKernel}