\name{kernel.lin}
\alias{kernel.lin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear Kernel
}
\description{
Creates linear kernel. 
}
\usage{
kernel.lin(GWASdata, pathway, parallel = c("none", "cpu", "gpu"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GWASdata}{
  Object of the class \code{GWASdata}.
}
  \item{pathway}{
  Object of the class \code{pathway}.
}
  \item{parallel}{
  Use \code{"none"}, \code{"cpu"} or \code{"gpu"}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns object of class kernel. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Stefanie Friedrichs with contributions by Saskia Freytag and Ngoc-Thuy Ha
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GWASdata}}, \code{\link{pathway}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (GWASdata, pathway, parallel = c("none", "cpu", "gpu")) 
{
    SNPset <- GWASdata@anno$rsNumber[which(GWASdata@anno$Pathway == 
        pathway@id)]
    z <- as(GWASdata@geno[, as.character(SNPset)], "matrix")
    if (parallel == "none") {
        k <- tcrossprod(z)
    }
    if (parallel == "cpu") {
        print("sorry, not yet defined")
    }
    if (parallel == "gpu") {
        z <- magma(z, gpu = TRUE)
        k <- tcrossprod(z)
    }
    return(kernel(type = "linear", kernel = k, GWASdata = GWASdata, 
        pathway = pathway))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
