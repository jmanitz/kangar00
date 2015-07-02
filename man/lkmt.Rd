\name{lkmt}
\alias{lkmt}

\title{
Logistic kernel maschine test functions to evaluate kernel similarity matrices in a variance component test. 
}
\description{

}
\usage{
lkmt(formula, kernel, GWASdata, ...)
}
\arguments{
  \item{formula}{
a formular giving the regression nullmodel for the test. Notation is similiar to lm() or glm() functions. 
}
  \item{kernel}{
kernel square matrix of dimension equal to the number individuals
}
  \item{GWASdata}{
  Object of the class \code{GWASdata} containing the genotype and annotation files.
}
}

\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner
}

\keyword{lkmt}

