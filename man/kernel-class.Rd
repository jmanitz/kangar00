\name{kernel-class}
\Rdversion{1.1}
\docType{class}
\alias{kernel-class}
\alias{plot,kernel-method}
\alias{show,kernel-method}
\alias{summary,kernel-method}

\title{Class \code{kernel}}
\description{
This class describes different kind of kernels. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{kernel(...)}.
}
\section{Slots}{
  \describe{
    \item{\code{type}:}{Use "lin" for linear kernels, "adj" for size-adjusted or "net" for network-based kernels.}
    \item{\code{kernel}:}{Kernel matrix of dimension equal to individuals.}
    \item{\code{pathway}:}{Object of class \code{pathway}.}
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(object = "kernel")}: ... }
    \item{show}{\code{signature(object = "kernel")}: ... }
    \item{summary}{\code{signature(object = "kernel")}: Gives the number of individuals and genes for pathway.}
	 }
}
\seealso{
\code{\link{GWASdata}}, \code{\link{pathway}}
}
\examples{
showClass("kernel")
}
\keyword{classes}
\keyword{kernelClass}

