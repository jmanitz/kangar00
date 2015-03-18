\name{kernel.adj}
\alias{kernel.adj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Size-adjusted kernel
}
\description{
Calculates the size-adjusted kernel.
}
\usage{
kernel.adj(GWASdata, pathway)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GWASdata}{
  Object of the class \code{GWASdata}. 
}
  \item{pathway}{
  Object of the class \code{pathway}
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{GWASdata}}, \code{\link{pathway}}, \code{\link{kernel}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (GWASdata, pathway) 
{
    dat <- GWASdata@geno
    anno <- GWASdata@anno
    n <- pathway@id
    genes <- anno[anno[, 1] == paste(n, sep = ""), c(2, 4)]
    genes <- as.character(genes[!duplicated(genes), 1])
    gene.counts <- table(genes)
    max.length <- max(gene.counts)
    kerneltimes <- matrix(rep(0, (dim(dat)[1] * dim(dat)[1])), 
        nrow = dim(dat)[1])
    effSNPs <- EffectiveNumberSNPs(GWASdata, genes.id = names(gene.counts[gene.counts >= 
        10]))
    max.eff.length <- max(as.numeric(effSNPs[, 2]), na.rm = T)
    for (g in names(gene.counts[gene.counts >= 10])) {
        z <- as(dat[, as.character(unique(anno[(anno[, 1] == 
            paste(n, sep = "")) & (anno[, 2] == g) & (!duplicated(anno[, 
            c(1, 2, 4)])), 4]))], "matrix")
        z <- z[, apply(z, 2, sum)/(2 * dim(z)[1]) >= 0.001]
        eff.length.gene <- as.numeric(effSNPs[which(effSNPs[, 
            1] == g), 2])
        distances <- matrix(rep(rowSums(z * z), dim(z)[1]), nrow = dim(z)[1]) - 
            2 * z \%*\% t(z) + matrix(rep(rowSums(z * z), dim(z)[1]), 
            nrow = dim(z)[1], byrow = TRUE)
        distances <- round(distances, digits = 3)
        length.gene <- as.numeric(dim(z)[2])
        alpha <- sqrt(eff.length.gene/max.eff.length)
        roh <- (mean(c(distances)))^(-alpha) * (eff.length.gene/length.gene)^(-alpha)
        kerneltimes <- kerneltimes + -roh * (distances/length.gene)^(alpha)
    }
    k <- exp(sqrt(1/(length(unique(genes)))) * kerneltimes)
    return(kernel(type = "size-adjusted", kernel = k, GWASdata = GWASdata, 
        pathway = pathway))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{kernel.adj}
\keyword{size adjused kernel}% __ONLY ONE__ keyword per line
