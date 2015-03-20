\name{definition}
\alias{definition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Number of genes of SNPs. 
}
\description{
Calculates with help of annotation file number of genes and SNPs. 
}
\usage{
definition(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Object of class \code{GWASdata}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{tab}{
  Returns vector which contains number of genes and SNPs.
  }
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
\code{\link{GWASdata}} for more information about \code{object}. 
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object) 
{
    anno <- object@anno
    anno <- anno[!duplicated(anno[, 1:4]), ]
    nrsnps <- table(anno[, 1])
    anno <- anno[!duplicated(anno[, 1:2]), ]
    nrgenes <- table(anno[, 1])
    tab <- cbind(nrgenes, nrsnps)
    return(tab)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{definition}
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
