\name{kernel.lin}
\alias{kernel.lin}
\title{
Linear Kernel
}
\description{
Creates linear kernel. 
}
\usage{
kernel.lin(data, pathway, parallel = c("none", "cpu", "gpu"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  Object of the class \code{GWASdata}.
}
  \item{pathway}{
  Object of the class \code{pathway}.
}
  \item{parallel}{
  Use \code{"none"}, \code{"cpu"} or \code{"gpu"}. 
}
  \item{...}{
  ???
  }
}

\value{
Returns object of class kernel. 
}

\author{
Stefanie Friedrichs with contributions by Saskia Freytag and Ngoc-Thuy Ha
}

\seealso{
\code{\link{GWASdata}}, \code{\link{pathway}}
}

\keyword{LinearKernel}