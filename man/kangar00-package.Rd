\name{kangar00-package}
\alias{kangar00-package}
\alias{kangar00}
\docType{package}
\title{
kangar00 package
}
\description{
This package includes methods to extract information on pathways, genes and SNPs from 
online databases and to evaluate these data using the logistic kernel machine test 
(LKMT) (Liu et al. 2008).  \cr
We defined SNP sets representing genes and whole pathways using knowledge on
 gene membership and interaction from the Kyoto Encyclopedia of Genes and 
 Genomes (KEGG) database (Kanehisa et al. 2014). SNPs are mapped to genes via 
  base pair positions of SNPs and transcript start and end points of genes as 
  documented in the Ensemble database (Cunningham et al. 2015).  \cr
In the LKMT, we employed the linear kernel (Wu et al. 2010) as well as two more 
advanced kernels,  adjusting for size bias in the number of SNPs and genes in a
 pathway and incorporating the network structure of genes within the pathway, 
 respectively (Freytag et al. 2012, 2014).  \cr
  P-values are derived in a variance component test using a moment matching methode (Schaid, 2010) or Davies' algorithm (Davies, 1980). 
}
\details{
\tabular{ll}{
Package: \tab kangar00\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-04\cr
License: \tab GPL-2 \cr
Depends: \tab methods, igraph, biomaRt, KEGGgraph, sqldf \cr
}
}
\author{
Juliane Manitz, Stefanie Friedrichs, Patricia Burger, Benjamin Hofner  \cr
Maintainer: Stefanie Friedrichs <stefanie.friedrichs@med.uni-goettingen.de>
}
\references{
\itemize{
 \item Schaid DJ: Genomic similarity and kernel methods I: advancements by building on mathematical and statistical foundations. Hum Hered 2010, 70:109-131.   \cr
 \item	Davies R: Algorithm as 155: the distribution of a linear combination of chi-2 random variables. J R Stat Soc Ser C 1980, 29:323-333.  \cr
 \item Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin X: Powerful SNP-Set Analysis for Case-Control Genome-Wide Association Studies. Am J Hum Genet 2010, 86:929-42
 \item	Fiona Cunningham, M. Ridwan Amode, Daniel Barrell et al. Ensembl 2015. Nucleic Acids Research 2015 43 Database issue:D662-D669  \cr
 \item Kanehisa, M., Goto, S., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.; Data, information, knowledge and principle: back to metabolism in KEGG. Nucleic Acids Res. 42, D199-D205 (2014). \cr
 \item Freytag S, Bickeboeller H, Amos CI, Kneib T, Schlather M: A Novel Kernel for Correcting Size Bias in the Logistic Kernel Machine Test with an Application to Rheumatoid Arthritis. Hum Hered. 2012, 74(2):97-108. \cr
 \item Freytag S, Manitz J, Schlather M, Kneib T, Amos CI, Risch A, Chang-Claude J, Heinrich J, Bickeboeller H: A network-based kernel machine test for the identification of risk pathways in genome-wide association studies. Hum Hered. 2013, 76(2):64-75. \cr
 \item Liu D, Ghosh D, Lin X. Estimation and testing for the effect of a genetic pathway on a disease outcome using logistic kernel machine regression via logistic mixed models. BMC Bioinformatics. 2008 9:292. \cr
}
}

\keyword{ package }
\seealso{

}
\examples{

}
