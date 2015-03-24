\name{lkmt-class}
\Rdversion{1.1}
\docType{class}
\alias{lkmt-class}
\alias{show,lkmt-method}
\alias{summary,lkmt-method}

\title{Class \code{"lkmt"}}
\description{
Not yet. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("lkmt", ...)}.
 
}
\section{Slots}{
  \describe{
    \item{\code{GWASdata}:}{Object of class \code{"GWASdata"}.}
    \item{\code{df}:}{Degrees of freedom}
    \item{\code{formula}:}{Object of class \code{"formula"}: Gives formula for general linear model.}
    \item{\code{kernel}:}{Object of class \code{"kernel"}.}

    \item{\code{p.value}:}{Object of class \code{"vector"} which contains pvalues.}
    \item{\code{statistic}:}{Test statistic}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "lkmt")}: ... }
    \item{summary}{\code{signature(object = "lkmt")}: ... }
	 }
}


\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner
}

\examples{
showClass("lkmt")
}
\keyword{classes}
