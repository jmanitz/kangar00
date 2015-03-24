\name{EffectiveNumberSNPs}
\alias{EffectiveNumberSNPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Effective number of SNPs
}
\description{
This function calculates the effective number of SNPs. 
}
\usage{
EffectiveNumberSNPs(object, genes.ids)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Object of the type \code{GWASdata}. 
}
  \item{genes.ids}{
  Gene IDs 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{effectivesnps}{Effective number of SNPs.}
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
}
\examples{
## ORIGINAL CODE!!!!
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, genes.ids) 
{
    dat <- object@geno
    anno <- object@anno
    effective <- function(x) {
        length(x) * (1 - (length(x) - 1) * var(x)/(length(x)^2))
    }
    effectivesnps <- matrix(NA, nrow = length(genes.ids), ncol = 2)
    effectivesnps[, 1] <- as.character(genes.ids)
    for (g in genes.ids) {
        z <- as(dat[, as.character(unique(anno[(anno[, 2] == 
            g), 4]))], "matrix")
        z <- z[, apply(z, 2, sum)/(2 * dim(z)[1]) >= 0.001]
        eigen.values <- eigen(cor(z), symmetric = TRUE, only.values = TRUE)$values
        effectivesnps[effectivesnps[, 1] == g, 2] <- effective(eigen.values)
    }
    return(effectivesnps)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{EffectiveNumberSNPs}
