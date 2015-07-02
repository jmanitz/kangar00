\name{kernel.adj}
\alias{kernel.adj}
\title{
Size-adjusted kernel
}
\description{
Calculates the size-adjusted kernel which takes into consideration the numbers of SNPs and genes in a pathway to correct for size bias.
}
\usage{
kernel.adj(GWASdata, pathway)
}
\arguments{
  \item{GWASdata}{
  Object of the class \code{GWASdata}, containing the genotypes of the individuals for which a kernel will be calculated. 
}
  \item{pathway}{
  Object of the class \code{pathway}, containing information on the pathway for which a kernel will be calculated. 
}
  \item{parallel}{
  Argument specifying whether to calculate kernels on several CPUs in parallel (\code{"cpu"}) or on ond CPU (\code{"none"}) or on a GPU (\code{"gpu"}). 
}
}

\value{
Returns an object of class kernel, including the similaritymatrix of the corresponding pathway for the considered individuals.   
}
\references{
Freytag S, Bickeböller H, Amos CI, Kneib T, Schlather M: A Novel Kernel for Correcting Size Bias in the Logistic Kernel Machine Test with an Application to Rheumatoid Arthritis. Hum Hered. 2012, 74(2): 97–108.
}
\seealso{
\code{\link{GWASdata}}, \code{\link{pathway}}, \code{\link{kernel}}
}
\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner with contributions by Saskia Freytag and Ngoc-Thuy Ha.
}

\keyword{kernel.adj}
\keyword{SizeAdjusedKernel}
