\name{score.test}
\alias{score.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Score Test Using Sattherthwaite Approximation 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
score.test(kernels, nullmodel, pd.check = TRUE, kernel.out = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{kernels}{
  List of kernels-matrices, or a kernel-matrix
}
  \item{nullmodel}{
  A GLM-object.
}
  \item{pd.check}{
  Boolean, whether to check for positive definite or not.
}
  \item{kernel.out}{
  Boolean, whether the kernel need to be omitted.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Schaid, D. J. (2002) \emph{Score Tests for Association between Traits and Haplotypes when Linkage Phase Is Ambiguous} The American Journal of Human Genetics
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
function (kernels, nullmodel, pd.check = TRUE, kernel.out = FALSE) 
{
    if (is.matrix(kernels)) {
        kernels <- list(pathway1 = kernels)
    }
    else if (!is.list(kernels)) {
        stop("kernels should be either a kernel-matrix or a list of kernel matrices!")
    }
    if (sum(is(nullmodel) \%in\% c("glm", "lm")) == 0) {
        stop("nullmodel should be a glm- or lm-object!")
    }
    if (is.null(nullmodel$x)) {
        stop("The glm-object should have a design-matrix x!")
    }
    nas <- nullmodel$na.action
    pathwaynames <- names(kernels)
    Y <- nullmodel$y
    X <- nullmodel$x
    mui <- nullmodel$fitted.values
    W <- Diagonal(length(mui), mui * (1 - mui))
    WX <- W \%*\% X
    P <- W - WX \%*\% solve(t(X) \%*\% WX, t(X) \%*\% W)
    all.mod <- list()
    for (k in pathwaynames) {
        K <- kernels[[k]]
        K <- as.matrix(K)
        if (!is.null(nas)) {
            K <- K[-nas, -nas]
        }
        TT <- 1/2 * t(Y - mui) \%*\% K \%*\% (Y - mui)
        PK <- P \%*\% K
        A <- W \%*\% PK \%*\% P \%*\% W
        VarT <- 1/2 * sum(diag((PK \%*\% PK))) + 1/4 * sum(diag(A)^2 * 
            mui * (1 - mui) * (1 + 6 * mui^2 - 6 * mui))
        ExpT <- 1/2 * sum(diag((PK)))
        a <- VarT/(2 * ExpT)
        d <- (2 * ExpT^2)/VarT
        p.value <- pchisq(TT/a, d, lower.tail = FALSE)
        mod <- list(p.value = p.value, statistic = TT/a, parameter = d, 
            method = "Dan Schaid Score Test")
        names(mod$statistic) <- "Chi Squared Test Statictic Value"
        names(mod$parameter) <- "Degree of Freedom"
        if (kernel.out == T) {
            mod$kernel <- K
        }
        class(mod) <- "htest"
        all.mod[[as.character(k)]] <- mod
    }
    return(all.mod)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{score.test}
\keyword{schaid}
