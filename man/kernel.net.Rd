\name{kernel.net}
\alias{kernel.net}
\title{
Network-based kernel
}
\description{
Calculates a network-based kernel to be evaluated in the logistic kernel machine test. Not only information on gene membership and gene/pathway size in number of SNPs is incorporated in the calculation, but also the interaction structure of genes in the pathway. 
}
\usage{
kernel.net(data, pathway)
}

\arguments{
  \item{data}{
  Object of the class \code{GWASdata} containing the genotypes of the individuals for which a kernel will be calculated. 
}
  \item{pathway}{
  Object of the class \code{pathway} containing information on the pathway for which a kernel will be calculated.
}
  \item{parallel}{
  Argument specifying whether to calculate kernels on several CPUs in parallel (\code{"cpu"}) or on ond CPU (\code{"none"}) or on a GPU (\code{"gpu"}). 
}
}

\value{
Returns object of class kernel, including the similaritymatrix of the corresponding pathway for the considered individuals. 
}
\references{Freytag S, Manitz J, Schlather M, Kneib T, Amos CI, Risch A, Chang-Claude J, Heinrich J, Bickeböller H: A network-based kernel machine test for the identification of risk pathways in genome-wide association studies. Hum Hered. 2013, 76(2):64-75. 
}
\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner with contributions by Saskia Freytag and Ngoc-Thuy Ha.
}

\keyword{NetworkKernel}

