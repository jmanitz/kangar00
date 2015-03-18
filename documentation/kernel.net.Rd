\name{kernel.net}
\alias{kernel.net}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Network-based kernel
}
\description{
Calculates a network-based kernel
}
\usage{
kernel.net(GWASdata, pathway)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GWASdata}{
  Object of the class \code{GWASdata}.
}
  \item{pathway}{
  Object of the class \code{pathway}.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
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
    SNPset <- GWASdata@anno$rsNumber[which(GWASdata@anno$Pathway == 
        pathway@id)]
    Z <- as(GWASdata@geno[, as.character(SNPset)], "matrix")
    N <- as.matrix(pathway@adj)
    N[N != 0] <- pathway@sign
    all.genes <- data.frame(colnames(N))
    anno.genes <- data.frame(unique(anno[anno[, 1] == n, 2]), 
        0)
    detect <- merge(all.genes, anno.genes, by.x = 1, by.y = 1, 
        all.x = T)
    remov <- as.character(detect[is.na(detect[, 2]), 1])
    if (length(remov) != 0) {
        for (g in remov) {
            z <- which(colnames(N) == g)
            vec <- rbind(N[z, ], seq(1:length(N[z, ])))
            vec <- vec[, vec[1, ] != 0]
            if (length(vec) > 0) {
                for (i in 1:(length(vec[1, ]) - 1)) {
                  for (j in (i + 1):length(vec[1, ])) {
                    if (N[vec[2, i], vec[2, j]] != 0) {
                      print("Edge will be removed!")
                    }
                    N[vec[2, i], vec[2, j]] <- N[vec[2, i], vec[2, 
                      j]] + vec[1, i] * vec[1, j]
                    N[vec[2, j], vec[2, i]] <- N[vec[2, j], vec[2, 
                      i]] + vec[1, i] * vec[1, j]
                  }
                }
            }
        }
        N <- N[-z, -z]
    }
    N[N > 1] <- 1
    N[N < (-1)] <- -1
    diag(N) <- 1
    lambda <- min(eigen(N, only.values = TRUE, symmetric = TRUE)$values)
    if (lambda < 0) {
        rho <- 1/(1 - lambda)
        N <- rho * N + (1 - rho) * diag(dim(N)[1])
    }
    A <- matrix(data = 0, nrow = nrow(Z), ncol = nrow(N))
    colnames(A) <- rownames(N)
    rownames(A) <- as.character(unique(anno[anno[, 1] == paste(n, 
        sep = ""), 4]))
    for (i in 1:ncol(A)) {
        snpsi <- as.character(unique(anno[anno[, 2] == paste(colnames(A)[i], 
            sep = ""), 4]))
        A[snpsi, i] <- 1
    }
    A.star <- A
    for (i in 1:length(A[1, ])) {
        A.star[, i] <- A[, i]/sum(A[, i])
    }
    K <- Z \%*\% A.star \%*\% N \%*\% t(A.star) \%*\% t(Z)
    return(kernel(type = "network", kernel = K, GWASdata = GWASdata, 
        pathway = pathway))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
