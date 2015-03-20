\name{lkmt-class}
\Rdversion{1.1}
\docType{class}
\alias{lkmt-class}
\alias{show,lkmt-method}
\alias{summary,lkmt-method}

\title{Class \code{"lkmt"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lkmt", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{formula}:}{Object of class \code{"formula"}: Gives formula for general linear model.}
    \item{\code{kernel}:}{Object of class \code{"kernel"}.}
    \item{\code{GWASdata}:}{Object of class \code{"GWASdata"}.}
    \item{\code{pval}:}{Object of class \code{"vector"} which contains pvalues.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "lkmt")}: ... }
    \item{summary}{\code{signature(object = "lkmt")}: ... }
	 }
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("lkmt")
}
\keyword{classes}
