#' Example phenotype file for 50 individuals.
#'
#' A dataset containing simulated example phenotypes for 50 individuals 
#' Row-names include the identifiers of 50 example individuals.
#'
#' @format A \code{data frame} with 50 rows and 3 variables:
#' \describe{
#'  \item{pheno}{includes the case-control status for each individual, coded as
#'    1(case) or 0 (conrol)}
#'  \item{sex}{includes gender information for the 50 individuals, coded as 1 
#'    (male) or 0 (female)}
#'  \item{age}{numerical value giving the persons age}
#' }
#' @source simulated data
"pheno"

#' Example annotation file for three pathways. 
#'
#' A dataset containing an annotation example for 4056 SNPs in three different 
#' pathways. 
#'
#' @format A \code{data frame} with 4056 rows and 5 variables:
#' \describe{
#'  \item{pathway}{includes KEGG identifiers of three example pathways}
#'  \item{gene}{names of genes in the pathways}
#'  \item{chr}{specifies the chromosom}
#'  \item{snp}{includes rs-numbers of example SNPs}
#'  \item{position}{gives positions of example SNPs}
#' }
#' @source simulated data
"anno"

#' Example genotypes for 50 individuals.
#'
#' A matrix containing example genotypes for 4056 SNPs of 50 individuals. 
#' Column-names give the rs-numbers of 4056 example SNPs, row-names the 
#' identifiers of 50 example individuals.
#'
#' @format A \code{matrix} with 5 rows and 4056 columns:
#' \describe{
#'  each entry in the matrix represents a simulated minor allele count for the
#' corresponding SNp and individual.
#' }
#' @source simulated data
"geno"

#' Example GWASdata object.
#'
#' An object of type GWASdata  dataset containing the example files for 
#' annotatio, phenotypes and genotypes. 
#'
#' @format An object of class \code{GWASdata}:
#' \describe{
#'  \item{geno}{contains example genotypes}
#'  \item{anno}{example annotation for threee pathways}
#'  \item{pheno}{exemplary genotypes for the iondividuals listed in geno}
#'  \item{desc}{a description of the GWAS study, here 'example study'}
#' }
#' @source simulated data
"gwas"

#' Example network-based kernel matrix for pathway hsa04020.
#'
#' An example of a kernel object.  
#'
#' @format An object of class \code{kernel} and type 'network' for the pathway 
#' hsa04020.
#' \describe{
#'  \item{type}{specifies which kernelfunction was used to calculate the kernel}
#' \item{kernel}{includes the kernelmatrix calculated for the pathway}
#' \item{pathway}{includes the pathway object of the pathway, for which the 
#' kernelmatrix was calculated}
#' }
#' @source simulated data and Ensembl extract
"net.kernel.hsa04020"

#' Example test result for the network-based kernel for pathway hsa04020.
#'
#' An object of class \code{lkmt} containing exemplary test results for an 
#' application of the logistic kernel machine test, derived from the example data.  
#'
#' @format An object of class \code{lkmt} for the network-based kernel and 
#' the pathway hsa04020.
#' \describe{
#'  \item{formular}{gives the formular defining the nullmodel used in the 
#' logistic kernel machine test}
#' \item{kernel}{includes the \code{kernel} object of the pathway to be evaluated}
#' \item{GWASdata}{gives the \code{GWASdata} object including the study data 
#' considered in testing}
#' \item{statistic}{gives the value of the test statistic}
#' \item{df}{specifies the degrees of freedom}
#' \item{p.value}{includes teh p-value resulting from the test}
#' }
#' @source simulated data and Ensembl extract
"lkmt.net.kernel.hsa04020"

#' Example \code{pathway} object for pathway hsa04020. 
#'
#' An object of class \code{pathway} for the pathway with KEGG identifier hsa04020.  
#'
#' @format A \code{pathway} object including 180 genes.
#' \describe{
#'  \item{id}{KEGG identifier of the example pathways}
#'  \item{adj}{gives the quadratic adjacency matrix for the pathway and with 
#' that the network topology. Matrix dimensionequal the number of genes in the 
#' pathway}
#'  \item{sign}{includes a vector of signs to distinguish activations and 
#' inhibitions in the adjacency matrix }
#' }
#' @source simulated data and Ensembl extract
"hsa04020"

#' Example \code{pathway_info} object for pathway hsa04022. 
#'
#' An object of class \code{pathway_info} for the pathway with KEGG 
#' identifier hsa04020.  
#'
#' @format A \code{pathway_info} object including information on 163 genes.
#' \describe{
#'  \item{info}{a \code{data frame} including information on genes included in 
#' pathway. Has columns 'pathway', 'gene_start', 'gene_end', 'chr', and 'gene'}
#' }
#' @source Ensembl extract
"hsa04022_info"

#' Example \code{snp_info} object for SNP rs10243170. 
#'
#' An object of class \code{snp_info} for rs10243170. 
#'
#' @format A \code{snp_info} object including information on the SNP as 
#' extracted from the Ensembl database.
#' \describe{
#'  \item{info}{a \code{data frame} including the extracted information on the 
#' SNP. Columns given are 'chr', 'position', and 'rsnumber'}
#' }
#' @source Ensembl extract
"rs10243170_info"
