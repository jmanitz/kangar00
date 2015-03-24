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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{score.test}
\keyword{schaid}
