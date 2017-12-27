#' Example phenotype file for 50 individuals.
#'
#' A dataset containing simulated example phenotypes for 50 individuals 
#' row names include the identifiers of 50 example individuals.
#'
#' @format A \code{data frame} with 50 rows and 3 variables:
#' \describe{
#'  \item{pheno}{includes the case-control status for each individual, coded as
#'    1(case) or 0 (control)}
#'  \item{sex}{includes gender information for the 50 individuals, coded as 1 
#'    (male) or 0 (female)}
#'  \item{age}{numerical value giving the persons age}
#' }
#' @usage data(pheno)
#' @examples
#' data(pheno)
#' head(pheno)
#' # create gwas object
#' data(geno)
#' data(anno)
#' gwas <- new('GWASdata', pheno=pheno, geno=geno, anno=anno, desc="some study") 
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
#'  \item{chr}{specifies the chromosome}
#'  \item{snp}{includes rs-numbers of example SNPs}
#'  \item{position}{gives positions of example SNPs}
#' }
#' @usage data(anno)
#' @examples
#' data(anno)
#' head(anno)
#' # create gwas object
#' data(pheno)
#' data(geno)
#' gwas <- new('GWASdata', pheno=pheno, geno=geno, anno=anno, desc="some study") 
#' @source simulated data
"anno"

#' Example genotypes for 50 individuals.
#'
#' A matrix containing example genotypes for 4056 SNPs of 50 individuals. 
#' Column names give the rs-numbers of 4056 example SNPs, row names the 
#' identifiers of 50 example individuals.
#'
#' @format A \code{matrix} with 5 rows and 4056 columns:
#' \describe{
#'  each entry in the matrix represents a simulated minor allele count for the
#' corresponding SNP and individual.
#' }
#' @usage data(geno)
#' @examples
#' data(geno)
#' head(geno)
#' # create gwas object
#' data(pheno)
#' data(anno)
#' gwas <- new('GWASdata', pheno=pheno, geno=geno, anno=anno, desc="some study") 
#' @source simulated data
"geno"

#' Example \code{\link{GWASdata}} object.
#'
#' An object of type GWASdata containing the example files for 
#' annotation, phenotypes and genotypes. 
#'
#' @format An object of class \code{\link{GWASdata}}:
#' \describe{
#'  \item{geno}{contains example genotypes}
#'  \item{anno}{example annotation for three pathways}
#'  \item{pheno}{exemplary phenotypes for all 'genotyped' individuals}
#'  \item{desc}{a description of the GWAS study, here 'example study'}
#' }
#' @usage data(gwas)
#' @examples
#' # create gwas object
#' data(pheno)
#' data(geno)
#' data(anno)
#' gwas <- new('GWASdata', pheno=pheno, geno=geno, anno=anno, desc="some study") 
#' @source simulated data
"gwas"

#' Example network-based kernel matrix for pathway hsa04020.
#'
#' An example of a kernel object.  
#'
#' @format An object of class \code{\link{kernel}} and type 'network' for the pathway 
#' hsa04020.
#' \describe{
#'  \item{type}{specifies which kernel function was used to calculate the kernel}
#' \item{kernel}{includes the kernel matrix calculated for the pathway}
#' \item{pathway}{includes the \code{\link{pathway}} object of the pathway, for which  
#' the kernel matrix was calculated}
#' }
#' @usage data(net.kernel.hsa04020)
#'
#' @examples
#' data(net.kernel.hsa04020)
#' # derivation 
#' data(gwas)
#' data(hsa04020)
#' net_kernel <- calc_kernel(gwas, hsa04020, knots=NULL, type='net', calculation='cpu')
#' # are the value differences smaller than machine epsilon?
#' all(abs(net.kernel.hsa04020@kernel - net_kernel@kernel) < sqrt(.Machine$double.eps))
#' 
#' @source simulated data and Ensembl extract
"net.kernel.hsa04020"

#' Example test result for the network-based \code{\link{kernel}} for 
#' \code{\link{pathway}} hsa04020.
#'
#' An object of class \code{\link{lkmt}} containing exemplary test results for an 
#' application of the logistic kernel machine test, derived from the example data.  
#'
#' @format An object of class \code{\link{lkmt}} for the network-based  
#' \code{\link{kernel}} and the \code{\link{pathway}} hsa04020.
#' \describe{
#'  \item{formular}{gives a formular defining the nullmodel used in the 
#' logistic kernel machine test}
#' \item{kernel}{includes the \code{\link{kernel}} object of the 
#' \code{\link{pathway}} to be evaluated}
#' \item{GWASdata}{gives the \code{\link{GWASdata}} object including the study data 
#' considered in testing}
#' \item{statistic}{gives the value of the test statistic}
#' \item{df}{specifies the degrees of freedom}
#' \item{p.value}{includes teh p-value resulting from the test}
#' }
#' @usage data(lkmt.net.kernel.hsa04020)
#' @source simulated data and Ensembl extract
"lkmt.net.kernel.hsa04020"

#' Example \code{\link{pathway}} object for pathway hsa04020. 
#'
#' An object of class \code{\link{pathway}} for the pathway with KEGG 
#' identifier hsa04020.  
#'
#' @format A \code{\link{pathway}} object including 180 genes.
#' \describe{
#'  \item{id}{KEGG identifier of the example pathways}
#'  \item{adj}{gives the quadratic adjacency \code{matrix} for the pathway and with 
#' that the network topology. Matrix dimensions equal the number of genes in the 
#' pathway}
#'  \item{sign}{includes a \code{vector} of signs to distinguish activations and 
#' inhibitions in the adjacency \code{matrix} }
#' }
#' @usage data(hsa04020)
#' @source simulated data and Ensembl extract
"hsa04020"

#' Example \code{\link{pathway_info}} object for \code{\link{pathway}} hsa04022. 
#'
#' An object of class \code{\link{pathway_info}} for the \code{\link{pathway}}
#' with KEGG identifier hsa04020.  
#'
#' @format A \code{\link{pathway_info}} object including information on 163 genes.
#' \describe{
#'  \item{info}{a \code{data frame} including information on genes included in 
#' pathway. Has columns 'pathway', 'gene_start', 'gene_end', 'chr', and 'gene'}
#' }
#' @usage data(hsa04022_info)
#' @examples
#' \dontrun{
#' pathway_info('hsa04020')
#' }
#' @source Ensembl extract
"hsa04022_info"

#' Example \code{\link{snp_info}} object for SNP rs10243170. 
#'
#' An object of class \code{\link{snp_info}} for rs10243170. 
#'
#' @format A \code{\link{snp_info}} object including information on the SNP as 
#' extracted from the Ensembl database.
#' \describe{
#'  \item{info}{a \code{data frame} including the extracted information on the 
#' SNP. Columns given are 'chr', 'position', and 'rsnumber'}
#' }
#' @usage data(rs10243170_info)
#' @examples
#' snp_info("rs10243170")
#'
#' @source Ensembl extract
"rs10243170_info"
