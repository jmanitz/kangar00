\name{kernel}
\alias{kernel}
\title{
New kernel
}
\description{
This function creates a new kernel. 
}
\usage{
kernel(type, kernel, GWASdata, pathway)
}
\arguments{
    \item{\code{type}:}{Use "lin" for linear kernels, "adj" for size-adjusted or "net" for network-based kernels.}
    \item{\code{kernel}:}{Kernel matrix of dimension equal to individuals.}
    \item{\code{GWASdata}:}{Object of class \code{GWASdata}.}
    \item{\code{pathway}:}{Object of class \code{pathway}.}
}

\seealso{
\code{\link{kernel-class}}
}
\examples{

}
\keyword{kernel}

