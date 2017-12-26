#################################################
#
# kernel object functions
#
#################################################

#' @include GWASdata.r 
#' @include pathway.r
NULL

#######  kernel class definition  ################
#' An S4 class representing a kernel matrix calculated for a pathway
#'
#' @rdname kernel-class
#'
#' @slot type A \code{character} representing the kernel type: Use 
#' \code{'lin'} for linear kernel, \code{'sia'} for the size-adjusted 
#' or \code{'net'} for the network-based kernel.
#' @slot kernel A kernel \code{matrix} of dimension equal to the number of individuals
#' @slot pathway A \code{\link{pathway}} object
#'
#' @examples
#' data(gwas)
#' data(hsa04020)
#' net_kernel <- calc_kernel(gwas, hsa04020, knots=NULL, type='net', calculation='cpu')
#'
#' @author Juliane Manitz
#' @export
#' @import methods
kernel <- setClass('kernel',
                   representation(type='character', kernel='matrix',
                                  pathway='pathway'))
setValidity('kernel', function(object){
    msg  <- NULL
    valid <- TRUE
#    ## if knots are specified, skip validy check
#    further_args <- list(...)
#    if (!is.null(further_args) && "knots" %in% names(further_args)) {
#  	return(TRUE)
#    }
    # kernel matrix quadratic
#    if(nrow(object@kernel)!=ncol(object@kernel)){
#        valid <- FALSE
#        msg <- c(msg, "kernel matrix has to be quadratic")
#    }
    #	# nrow/ncol = # individuals
    #	if(nrow(object@kernel)!=nrow(object@GWASdata@pheno)){
    #	  valid <- FALSE
    #     msg <- c(msg, "kernel matrix has to be dimension equal to individuals in GWAS")
    #        }
    # isSymmetric(kernel)
    if( !isSymmetric(round(object@kernel),10) ){
        valid <- FALSE
        msg <- c(msg, "kernel matrix has to be symmetric")
    }
    # pos semi def?
    lambda <- min(eigen(object@kernel, only.values=TRUE, symmetric=TRUE)$values)
    if(lambda < -1e-6){ # smallest eigenvalue negative = not semipositive definite
        valid <- FALSE
        msg <- c(msg, "kernel matrix has to be positive semi-definite")
    }
    if(valid) TRUE else msg
})

#' An S4 class to represent a low-rank kernel for a SNPset at specified knots
#'
#' @rdname lowrank_kernel-class
#'
#' @slot type character, kernel type: Use \code{'lin'} for the linear kernel,
#' \code{'sia'} for the size-adjusted or \code{'net'} for the network-based kernel.
#' @slot kernel kernel \code{matrix} of dimension equal to individuals
#' @slot pathway \code{\link{pathway}} object
#'
#' @details This kernel is used for predictions. If observations and knots are 
#' equal, better construct a full-rank kernel of class \code{\link{kernel}}.
#'
#' @examples
#' data(gwas)
#' calc_kernel(gwas, hsa04020, knots=gwas, type='lin', calculation='cpu') 
#' \dontrun{
#' gwas2 <- new('GWASdata', pheno=pheno[1:10,], geno=geno[1:10,], anno=anno, desc=" study 2")
#' calc_kernel(gwas, hsa04020, knots = gwas2, type='net', calculation='cpu')
#' }
#' @author Juliane Manitz
#' @export
#' @import methods
lowrank_kernel <- setClass('lowrank_kernel',
                   slots=c(type='character', kernel='matrix', pathway='pathway'))

setValidity('lowrank_kernel', function(object){
    msg   <- NULL
    valid <- TRUE
    if( sum(is.na(object@kernel))>0 ){
        valid <- FALSE
        msg   <- c(msg, "kernel matrix has missing values")
    }
    if(valid) TRUE else msg
})
#################################### kernel object constructor #################

# calculate kernel object
setGeneric('calc_kernel', function(object, ...) standardGeneric('calc_kernel'))
#' Calculates the kernel-matrix for a pathway
#'
#' Uses individuals' genotypes to create a \code{\link{kernel}} object including 
#' the calculated kernel \code{matrix} for a specific \code{\link{pathway}}.
#' Each numeric value within this \code{matrix} is calculated
#' from two individuals' genotypevectors of the SNPs within 
#' the \code{\link{pathway}} by a kernel function. It can be interpreted as the genetic 
#' similiarity of the individuals. Association between the \code{\link{pathway}} and a 
#' binary phenotype (case-control status) can be evaluated
#' in the logistic kernel machine test, based on the \code{\link{kernel}} object. 
#' Three kernel functions are available. 
#'
#' @param object \code{GWASdata} object containing the genotypes of the 
#' individuals for which a \code{\link{kernel}} will be calculated.
#' @param pathway object of the class \code{\link{pathway}} specifying the SNP set 
#' for which a \code{\link{kernel}} will be calculated.
#' @param type \code{character} indicating the \code{\link{kernel}} type: 
#' Use \code{'lin'} to specify the linear kernel, \code{'sia'} for the 
#' size-adjusted or \code{'net'} for the network-based kernel.
#' @param knots \code{GWASdata} object, if specified a \code{\link{kernel}}  
#' will be computed.
#' @param calculation \code{character} specifying if the kernel \code{matrix} 
#' is computed on CPU or GPU. 
#' @param ... further arguments to be passed to \code{\link{kernel}} computations.
#'
#' @return Returns an object of class \code{\link{kernel}}, including the similarity 
#' \code{matrix} of the \code{\link{pathway}} for the considered individuals. \cr
#' If \code{knots} are specified low-rank kernel of class a \code{lowrank_kernel} 
#' will be returned, which is not necessarily quadratic and symmetric.
#' @details
#' Different types of kernels can be constructed:
#' \itemize{
#'   \item \code{type='lin'} creates the linear kernel assuming additive SNP 
#'effects to be evaluated in the logistic kernel machine test.
#'   \item \code{type='sia'} calculates the size-adjusted kernel which takes 
#' into consideration the numbers of SNPs and genes in a \code{\link{pathway}} 
#' to correct for size bias.
#'   \item \code{type='net'} calculates the network-based kernel. Here not only information on gene membership and gene/pathway size in number of SNPs is incorporated, but also the interaction structure of genes in the \code{\link{pathway}}.
#' }
#' For more details, check the references.
#' @references
#' \itemize{
#'  \item Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin X Powerful SNP-Set Analysis for Case-Control Genome-Wide Association Studies. Am J Hum Genet 2010, 86:929-42
#'  \item Freytag S, Bickeboeller H, Amos CI, Kneib T, Schlather M: A Novel Kernel for Correcting Size Bias in the Logistic Kernel Machine Test with an Application to Rheumatoid Arthritis. Hum Hered. 2012, 74(2):97-108.
#'  \item Freytag S, Manitz J, Schlather M, Kneib T, Amos CI, Risch A, Chang-Claude J, Heinrich J, Bickeboeller H: A network-based kernel machine test for the identification of risk pathways in genome-wide association studies. Hum Hered. 2013, 76(2):64-75.
#' }
#'
#' @examples
#' data(gwas)
#' data(hsa04020)
#' data(gwas)
#' net_kernel <- calc_kernel(gwas, hsa04020, knots=NULL, type='net', calculation='cpu')
#' show(net_kernel)
#' summary(net_kernel)
#' plot(net_kernel, hclust=TRUE)
#' 
#' @author Stefanie Friedrichs, Juliane Manitz
#' @rdname calc_kernel
#' @name calc_kernel
#' @aliases calc_kernel,GWASdata-method
#' @export
#' @seealso \code{\link{kernel-class}},\code{\link{pathway}}
setMethod('calc_kernel', signature(object = 'GWASdata'),
       definition = function(object, pathway, knots = NULL,
			                       type = c('lin', 'sia', 'net'),
                             calculation = c('cpu', 'gpu'), ...) {
           # user inputs
           type     <- match.arg(type)
           calculation <- match.arg(calculation)
	   if(!inherits(object, "GWASdata")){
	       stop("GWASdata must inherit from class 'GWASdata'")
	   }
	   if(!inherits(pathway, "pathway")){
	       stop("GWASdata must inherit from class 'pathway'")
	   }
	   if(!is.null(knots) && !inherits(knots, "GWASdata")){
	       stop("knots must inherit from class 'GWASdata'")
	   }
     
	   pwIdTest <- try(unique(object@anno$snp[which(object@anno$pathway == pathway@id)]))
     if(class(pwIdTest) == "try-error"){
       print(setdiff(pathway@id, object@anno$pathway))
       stop("The above-mentioned pathways cannot be found in the annotation file. \n")
     }
     
	   # transfer to specific kernel function
           k <- eval(parse(text=paste(type, "_kernel(
                           object = object, pathway = pathway,
                           knots = knots, calculation = calculation, ...)",sep='')))
           return(k)
})

############################### kernel functions ##############################
# calculate linear kernel
setGeneric('lin_kernel', function(object, ...) standardGeneric('lin_kernel'))
#' @describeIn calc_kernel Calculates a linear kernel 
# @export
setMethod('lin_kernel', signature(object = 'GWASdata'),
          definition = function(object, pathway, knots=NULL,
                       calculation = c('cpu', 'gpu'), ...) {
    calculation <- match.arg(calculation)
    lowrank <- !is.null(knots)
#    further_args <- list(...)
#    if (!is.null(further_args))
#        stop("handling of '...' not yet implemented")
    ## which SNPs are in specified pathway
    SNPset <- unique(object@anno$snp[which(object@anno$pathway == pathway@id)])
    #subset genotype data for specified SNP set
    Z1 <- as(object@geno[,as.character(SNPset)],'matrix')
    if(any(is.na(Z1)))
        stop("genotype information contains missing values")
    if(lowrank){
        Z2 <- knots@geno
        Z2 <- as(Z2[,as.character(SNPset)],'matrix')
        k <- Z1 %*% t(Z2)
        return(new('lowrank_kernel', type='lin', kernel=k, pathway=pathway))
    }else{
      # K=ZZ' kernel matrix = genetic similarity
      if(calculation=='cpu'){
        k <- tcrossprod(Z1)
      }
      if(calculation=='gpu'){   # Suggests gputools #
      if(require(gputools)){
        Z <- as.numeric(Z1)
        k <- gputools::gpuMatMult(Z1,t(Z1))
      }else{
        stop("Please install package 'gputools' to run matrix multiplication on GPU")
        }
      }
      #return kernel object
      return(new('kernel', type='lin', kernel=k, pathway=pathway))
    }
})


# create size-adjusted kernel
setGeneric('sia_kernel', function(object, ...) standardGeneric('sia_kernel'))
#' @describeIn calc_kernel Calculates a size adjusted-kernel
# @export 
setMethod('sia_kernel', signature(object = 'GWASdata'),
          definition = function(object, pathway, knots=NULL,
                       calculation = c('cpu', 'gpu'), ...) {
   
   calculation <- match.arg(calculation)
   if(calculation=='gpu'){ 
      stop("GPU calculation is not available for SIA kernel") 
   }
     
   # get Z1 matrix (old observations=full genotypes):    
   SNPset <- unique(object@anno$snp[which(object@anno$pathway==pathway@id)])
   z1 <- as(object@geno[,as.character(SNPset)],'matrix')
   if(any(is.na(z1)))
      stop("genotype information contains missing values")  
     
    #check if knots are specified
    lowrank <- !is.null(knots)  
      
    if(lowrank){
    #get Z2 matrix (new observations):
     z2 <- knots@geno
     z2 <- as(z2[,as.character(SNPset)],'matrix')
     if(any(is.na(z2)))
        stop("genotype information for new observations contains missing values")  
          
     anno  <- object@anno[object@anno[,"pathway"]==pathway@id, c("gene","snp")]  
     genes <- as.character(anno[,"gene"])
     gene.counts <- table(genes)
     max.length  <- max(gene.counts)
      
     EffectiveNumberSNPs <- function(gene.ids){
       effectivesnps <- cbind(as.character(gene.ids),rep(NA, length(gene.ids)))
       for(g in gene.ids){
        snps <- as.matrix(anno[anno$gene==g,"snp"])  
        #eigen.value calculation for effective number of snps on Z1:
        ev <- eigen(stats::cor(z1), symmetric=T)$values        
        effectivesnps[effectivesnps[,1]==g,2] <- 
          length(ev)*(1-(length(ev)-1)*stats::var(ev)/(length(ev)^2))
       }
     return(effectivesnps)    
     }
     gene.ids <- names(gene.counts[gene.counts >= 2]) 
     effSNPs <- EffectiveNumberSNPs(gene.ids) 
     max.eff.length <- max(as.numeric(effSNPs[,2]), na.rm=T)
    
     kerneltimes <- matrix(rep(0,nrow(z2)*nrow(z1)),nrow=nrow(z2))
    
     for(g in gene.ids){    
       g.z1 <- as(z1[,unique(anno[anno[,"gene"]==g,"snp"])],"matrix") #all inds 
       g.z2 <- as(z2[,unique(anno[anno[,"gene"]==g,"snp"])],"matrix") #knots         
       n1 <- nrow(g.z1) 
       n2 <- nrow(g.z2) 
         
       distances <- matrix(rep(rowSums(g.z2^2),n1), nrow=n2) - 2*g.z2%*%t(g.z1) 
                    + matrix(rep(rowSums(g.z1^2),n2), nrow=n2, byrow=TRUE) #ok.
       distances <- round(distances, digits=3)  #n* x n
                 
       length.gene <- as.numeric(ncol(g.z1))#snps in gene
       eff.length.g <- as.numeric(effSNPs[which(effSNPs[,1]==g),2])
       delta <- sqrt(eff.length.g/max.eff.length)
       roh  <- mean(c(distances))^(-delta)*(eff.length.g/length.gene)^(-delta)
       kerneltimes <- kerneltimes - roh*(distances/length.gene)^(delta)
       }                
    k <- exp( sqrt(1/(length(unique(anno[,"gene"])))) * kerneltimes )           
    return(lowrank_kernel(type='size-adjusted', kernel=k, pathway=pathway))
    }else{  
    genemat <- function(g, object){
        z <- z1
        z <- z[, apply(z,2,sum)/(2*nrow(z)) >= 0.001 &  apply(z,2,sum)/(2*nrow(z)) < 1] 
        #only snps maf >= 0.1%
        e.val <- eigen(stats::cor(z), symmetric=TRUE, only.values=TRUE)$values
        nn    <- length(e.val)
        a <- matrix( rep(rowSums(z*z),nrow(z)),nrow=nrow(z))
        distances <- a -2*tcrossprod(z) + t(a)
        distances <- round(distances, digits=3)
        return( list(distances, (nn*(1-(nn-1)*stats::var(e.val)/(nn^2))), ncol(z)) )
    }
        #matrix, eff.length.gene, length.gene
    genemat2 <- function(l, max.eff){
        delta <- sqrt(l[[2]]/max.eff) #delta <- sqrt(eff.length.gene/max.eff)
        roh   <- (mean(c(l[[1]])))^(-delta)*(l[[2]]/l[[3]])^(-delta)
        return(-roh*(l[[1]]/l[[3]])^(delta))
    }

    anno <- object@anno[object@anno[,"pathway"]==pathway@id, c("gene","snp")] 
    gene.counts <- table(anno[,"gene"]) #counts number of different SNPs per gene
    g.10 <- names(gene.counts[gene.counts >= 2]) #genes with >= 2 snps

    #[[1]]:genematrix, [[2]]:eff.length.gene, [[3]]:length.gene
    liste <- lapply(g.10, genemat, object)
    get2    <- function(l){ return(l[[2]]) }
    max.eff <- max( unlist( lapply(liste, get2) ) )

    kerneltimes <- matrix( rep(0,(dim(object@pheno)[1])^2), nrow=dim(object@pheno)[1])
    kerneltimes <- Reduce('+', lapply(liste,genemat2,max.eff))
    k <- exp( sqrt(1/(length(unique(anno[,"gene"])))) * kerneltimes )              
   return(kernel(type='size-adjusted',kernel=k,pathway=pathway)) 
   }   
})


# calculate network-based kernel
setGeneric('net_kernel', function(object, ...) standardGeneric('net_kernel'))
#' @describeIn calc_kernel Calculates a network-based kernel 
# @export
setMethod('net_kernel', signature(object = 'GWASdata'),
          definition = function(object, pathway, knots=NULL,
                       calculation = c('cpu', 'gpu'), ...) {
    calculation <- match.arg(calculation)                   
    ## check if knots are specified
    lowrank <- !is.null(knots)

    #genotype matrix Z, which SNPs are in specified pathway
    SNPset <- unique(object@anno$snp[which(object@anno$pathway==pathway@id)])
    #subset genotype data for specified SNP set
    Z1 <- as(object@geno[,as.character(SNPset)],'matrix')
    if (any(is.na(Z1)))
        stop("genotype information contains missing values")
    # compute kernel
    ANA <- get_ana(object@anno, SNPset, pathway)
    ## if knots are specified
    if(lowrank){
        Z2 <- knots@geno
        Z2 <- as(Z2[,as.character(SNPset)],'matrix')
        K <- Z1 %*% ANA %*% t(Z2)
        return(lowrank_kernel(type='network', kernel=K, pathway=pathway))
    }
   if(calculation=='cpu'){    
     K <- Z1 %*% ANA %*% t(Z1)
   }  
   if(calculation=='gpu'){   # Suggests gputools #
     K1 <- gputools::gpuMatMult(Z1,ANA)
     K  <- gputools::gpuMatMult(K1,t(Z1))
   }  
   #return kernel object
   return(kernel(type='network',kernel=K,pathway=pathway))
})

################################## helper function #############################

setGeneric('rewire_network', function(object, ...) standardGeneric('rewire_network'))
#' Rewires interactions in a \code{\link{pathway}}, which go through a gene not 
#' represented by any SNPs in the considered \code{\link{GWASdata}} dataset.   
#'
#' @export
#' @author Stefanie Friedrichs, Juliane Manitz
#'
#' @param object \code{\link{pathway}} object which's network \code{matrix} will be rewired
#' @param x A \code{vector} of gene names, indicating which genes are not represented
#' by SNPs in the considered \code{\link{GWASdata}} object and will be removed  
#' @return A \code{\link{pathway}} object including the rewired network \code{matrix}
#'
#' @name rewire_network
#' @rdname rewire_network
#' @aliases rewire_network,pathway-method
#' @examples
#' data(hsa04020)
#' summary(hsa04020)
#' hsa04020_rewired <- rewire_network(hsa04020, x=c('ADCY3', 'CALML3','GNAQ'))
#' summary(hsa04020_rewired)
#' 
## @references TODO Newman?   
setMethod('rewire_network', signature(object = 'pathway'),
          definition = function(object, x) {
          
    remov <- x
    #exit if no genes have to be removed 
    if(length(remov)==0){ return(object) }          
          
    N <- as.matrix(object@adj)
    N[N!=0] <- object@sign
    if(any(is.na(N)))
        stop("network information contains missing values")      
        
    # identify genes that need to be carried forward to the subnetwork
    a <- (N[remov,]!=0)
    # can be vector or matrix -> make matrix to apply colsums
    if(is.null(dim(a))){
       a <- rbind(rep(0,length(a)),a)
    }
    ind_sub <- which(colSums(a)!= 0)
    # extract the subnetwork
    Nsub <- N[ind_sub, ind_sub]

    #if gene to remove had no connections
    if(is.null(dim(Nsub))){
       N <- N[!(colnames(N)%in%colnames(N)[remov]),!(colnames(N)%in%colnames(N)[remov])]
       object@adj  <- abs(N) 
       object@sign <- as.vector(N[N!=0])
       return(object)
    }

    # exclude self-interaction
    diag(Nsub) <- 0
    # if nullmatrix
    if(sum(Nsub!=0)==0){   
       N <- N[!(colnames(N)%in%colnames(N)[remov]),!(colnames(N)%in%colnames(N)[remov])]
       object@adj  <- abs(N) 
       object@sign <- as.vector(N[N!=0])
       return(object)
    }

    # calculate the two-step network
    Nsub2step <- Nsub %*% Nsub
    Nsub2step[Nsub2step>0] <-  1
    Nsub2step[Nsub2step<0] <- -1

    # check whether interaction types contradict
    check_contradicts <- ((Nsub != 0) & (Nsub + Nsub2step == 0))
    if(any(check_contradicts)){
       Nsub2step[(Nsub != 0) & (Nsub + Nsub2step == 0)] <- 0
       message('Interaction types contradict after rewiring: Edges removed.')
    }
    # replace the subnetwork in the adjacency matrix
    N[ind_sub,ind_sub] <- Nsub2step
    # remove the genes and return network
    
    N <- N[!(colnames(N)%in%colnames(N)[remov]),!(rownames(N)%in%rownames(N)[remov])]
    object@adj  <- abs(N) 
    object@sign <- as.vector(N[N!=0])
    return(object)
})

setGeneric('get_ana', function(x, ...) standardGeneric('get_ana'))
#' Produce middle part of network kernel (for internal use)
#'
## @export
#' @author Juliane Manitz, Saskia Freytag, Stefanie Friedrichs
#'
#' @param x \code{data.frame} with annotation information as returned from
#' \code{get_anno}
#' @param SNPset \code{vector} with SNPs to be analyzed
#' @param pathway \code{pathway} object
#' @return matrix ANA' for inner part of network kernel
#' @keywords internal
setMethod('get_ana', signature = 'data.frame',
          definition = function(x, SNPset, pathway){

    ### remove genes that have no SNPs (not in anno):
    # genes in pathway
    net_genes <- get_genes(pathway)
    # genes in annotation -> genes with SNPs in GWASdata
    anno_sub <- x[x$pathway==pathway@id,]
    anno_genes <- unique(anno_sub$gene)
    # pathway genes that are not in annotation
    remov <- which(! net_genes %in% anno_genes)
    # rewire network -> separate function
    
    pathway <- rewire_network(pathway, remov)
    N <- as.matrix(pathway@adj)
    N[N!=0] <- pathway@sign
    if(any(is.na(N)))
        stop("network information contains missing values")
        
    # include selfinteractions for main effects
    diag(N) <- 1
    ## make N positive semidefinit
    N <- make_psd(N)

    #A: SNP to gene mapping
    Atab <- table(anno_sub[c('snp','gene')])
    Amat <- as(Atab, 'matrix')
    A <- Amat[SNPset,rownames(N)]    #A is colnames(Z) x rownames(N)

    #A*: size-adjustement for no of SNPs in each gene
    A.star <- t(t(A) / sqrt(colSums(A)))

    return(A.star %*% N %*% t(A.star))
})

setGeneric('make_psd', function(x, ...) standardGeneric('make_psd'))
#' Adjust network \code{matrix} to be positive semi-definite
#'
#' @param x A \code{matrix} specifying the network adjacency matrix.
#' @param eps A \code{numeric} value, setting the tolance for smallest 
#' eigenvalue adjustment
#' @return The \code{matrix} \code{x}, if it is positive definite and the closest
#' positive semi-definite \code{matrix} if \code{x} is not positive semi-definite.
#'
#' @details For a \code{matrix} N, the closest positive semi-definite \code{matrix} is 
#' calculated as N* = rho*N + (1+rho)*I, where I is the identity \code{matrix}
#' and rho = 1/(1 - lambda) with lambda the smallest eigenvalue of N. 
#' For more details check the references.
#'
#' @references
#' \itemize{
#'  \item Freytag S, Manitz J, Schlather M, Kneib T, Amos CI, Risch A, Chang-Claude J, Heinrich J, Bickeboeller H: A network-based kernel machine test for the identification of risk pathways in genome-wide association studies. Hum Hered. 2013, 76(2):64-75.
#' }
#'
#' @export
## @rdname make_psd
## @name make_psd
#' @aliases make_psd matrix
#' @author Juliane Manitz, Saskia Freytag, Stefanie Friedrichs
setMethod('make_psd', signature = 'matrix',
          definition = function(x, eps = sqrt(.Machine$double.eps)) {

    lambda <- min(eigen(x, only.values = TRUE, symmetric = TRUE)$values)
    ## use some additional tolerance to ensure semipositive definite matrices
    lambda <- lambda - sqrt(.Machine$double.eps)
    # smallest eigenvalue negative = not semipositive definite
    if (lambda < -1e-10) {
        rho <- 1/(1-lambda)
        x <- rho * x + (1-rho) * diag(dim(x)[1])
        ## now check if it is really positive definite by recursively calling
        x <- make_psd(x)
    }
    return(x)
})

#################################### basic methods for kernel #################
# show method
#' \code{show} displays the \code{kernel} object briefly
#' @param object An object of class \code{kernel}
#'
#' @examples
#' data(gwas)
#' data(hsa04020)
#' data(gwas)
#' net_kernel <- calc_kernel(gwas, hsa04020, knots=NULL, type='net', calculation='cpu')
#' show(net_kernel)
#'
#' @export
#' @rdname kernel-class
setMethod('show', signature('kernel'),
          definition = function(object){
              cat('An object of class ', class(object), ' of type ', object@type,' for pathway ', object@pathway@id, '\n\n',sep='')
              print(object@kernel)
              invisible(NULL)
          })

# summary method
setGeneric('summary', function(object) standardGeneric('summary'))

#' \code{summary} generates a \code{kernel} object summary including the number of 
#' individuals and genes for the \code{\link{pathway}}
#' @examples
#' summary(net_kernel)
#'
#' @export
#' @rdname kernel-class
#' @aliases summary,kernel,ANY-method
setMethod('summary', signature='kernel',
          definition = function(object){
              cat('An object of class ', class(object), ' of type ', object@type,' for pathway ', object@pathway@id, ' with values: \n\n',sep='')
              cat(paste('\nNumber of Individuals:',dim(object@kernel)[1]),'\n')
              cat(paste('The pathway contains',dim(object@pathway@adj)[1],'genes.\n'))#,'genes and', SNPno[2], 'SNPs.\n'))
              invisible(NULL)
          })

# plot method
if (!isGeneric("plot")) setGeneric('plot')
#' \code{plot} creates an image plot of a \code{kernel} object
#'
#' @param x the \code{kernel} object to be plotted.
#' @param y missing (placeholder).
#' @param hclust \code{logical}, indicating whether a dendrogram should be added.
#' @param ... further arguments to be passed to the function.
#'
#' @import lattice
#' @export
#' @rdname kernel-class
#' @examples
#' summary(net_kernel)
#' plot(net_kernel, hclust=TRUE)
#'
#' @aliases plot,kernel,ANY-method
setMethod('plot', signature(x='kernel',y='missing'),
          function(x, y=NA, hclust=FALSE, ...){
              if(hclust) {                                    
                  stats::heatmap(x@kernel, symm=TRUE, col=rev(grDevices::heat.colors(n=20)), Colv=NA,labRow=NA,labCol=NA, main=list(paste('Genetic Similarity Kernel Matrix for Pathway',x@pathway@id), cex=1.4), ...)
              }else{                                                 
                  print(lattice::levelplot(x@kernel, col.regions=rev(grDevices::heat.colors(n=20)), drop.unused.levels=FALSE, scales=list(alternating=0), main=paste('Genetic Similarity Kernel Matrix for Pathway',x@pathway@id)), ...)
              }
              invisible(NULL)
          })
