\name{kernel}
\alias{kernel}
\title{
New kernel constructer.
}
\description{
This function creates a new kernel. 
}
\usage{
kernel(type, kernel, pathway)
}
\arguments{
  \item{type}{Specifies if kernel symilarity-matrix was calculated using the linear, size-adjusted or network kwenel function}
  \item{kernel}{kernel matrix of dimension equal to individuals}
  \item{GWASdata}{GWASdata object}
  \item{pathway}{pathway object, giving information on the pathway on which the kernel matrix was calculated}
}

\seealso{
\code{\link{kernel-class}}
}

\keyword{kernel}

