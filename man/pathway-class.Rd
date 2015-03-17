\name{pathway-class}
\Rdversion{1.1}
\docType{class}
\alias{pathway-class}
\alias{pathway}
\alias{plot,pathway-method}
\alias{show,pathway-method}
\alias{summary,pathway-method}

\title{Class \code{pathway}}
\description{
This class describes pathways.
}
\section{Objects from the Class}{
New objects can be created by calls of the form \code{pathway(...)}.
}
\section{Slots}{
  \describe{
    \item{\code{id}:}{Object of class \code{"character"} describing pathway id, e.g. hsa00100 }
    \item{\code{adj}:}{Adjacency matrix}
    \item{\code{sign}:}{Object of class \code{"vector"} containing interaction type for each link }
  }
}
\section{Methods}{
  \describe{
    \item{plot}{\code{signature(object = "pathway")}: ... }
    \item{show}{\code{signature(object = "pathway")}: Summarize pathway information and prints  pathway adjacency matrix.}
    \item{summary}{\code{signature(object = "pathway")}: Information about nodes, links, activations and inhibitions.}
	 }
}
\seealso{
\code{\link{pathway}}
}
\examples{
showClass("pathway")
}
\keyword{classes}
\keyword{pathway class}

