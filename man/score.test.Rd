\name{score.test}
\alias{score.test}
\title{
Score Test Using Sattherthwaite Approximation 
}
\description{
Executes Score Test Using Sattherthwaite Approximation by Schaid
}
\usage{
score.test(kernels, nullmodel, pd.check = TRUE, kernel.out = FALSE)
}

\arguments{
  \item{kernels}{
  kernel-matrix for the pathway to be tested. 
}
  \item{nullmodel}{
  A GLM-object stating the nullmodel without genetic influence.
}
  \item{pd.check}{
  Boolean, whether to check for positive definite or not.
}
  \item{kernel.out}{
  Boolean, whether the kernel need to be omitted.
}
}
\references{
Schaid, D. J. (2002) \emph{Score Tests for Association between Traits and Haplotypes when Linkage Phase Is Ambiguous} The American Journal of Human Genetics
}
\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner
}
\keyword{score.test}
\keyword{schaid}
