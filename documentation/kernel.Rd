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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
structure(function (...) 
new("kernel", ...), className = structure("kernel", package = "kangar00"), package = "kangar00", class = structure("classGeneratorFunction", package = "methods"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{kernel}

