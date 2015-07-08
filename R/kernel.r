#################################################
#
# kernel object functions
#
#################################################

#' An S4 class to represent a kernel for a SNPset
#'
#' @rdname kernel-class
#'
#' @slot type character, kernel type: Use \code{"lin"} for linear kernels, 
#' \code{"sia"} for size-adjusted or \code{"net"} for network-based kernels.
#' @slot kernel kernel matrix of dimension equal to individuals
#' @slot GWASdata GWASdata object
#' @slot pathway pathway object
#' 
#' @author Juliane Manitz
#' @export
kernel <- setClass('kernel',
                   slots=c(type='character', kernel='matrix', pathway='pathway'))

setValidity('kernel', function(object){
    msg  <- NULL
    valid <- TRUE
    # kernel matrix quadratic
    if(nrow(object@kernel)!=ncol(object@kernel)){
        valid <- FALSE
        msg <- c(msg, "kernel matrix has to be quadratic")
    }
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
    # pos def?
    lambda <- min(eigen(object@kernel, only.values=TRUE, symmetric=TRUE)$values)
    if(lambda< -1e-6){ # smallest eigenvalue negative = not semipositive definite
        valid <- FALSE
        msg <- c(msg, "kernel matrix has to be positive semi-definite")
    }
    if(valid) TRUE else msg
})

# show method
#' \code{show} displays the kernel object briefly
#' @param object kernel object
#' @export
#' @rdname kernel-class
#' @aliases show,kernel,ANY-method
setMethod('show', signature='kernel',
          definition = function(object){
              cat('An object of class ', class(object), ' of type ', object@type,' for pathway ', object@pathway@id, '\n\n',sep='')
              print(object@kernel)
              invisible(NULL)
          })

# summary method
setGeneric('summary', function(object, ...) standardGeneric('summary'))

#' \code{summary} generates a kernel object summary including the number of individuals and genes for the pathway
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

#' \code{plot} creates an image plot of a kernel object
#'
#' @param y missing (placeholder)
#' @param hclust \code{logical}, indicating whether a dendrogram should be added
#' 
#' @export
#' @rdname kernel-class 
#' @aliases plot,kernel,ANY-method
setMethod('plot', signature(x='kernel',y='missing'),
          function(x, y=NA, hclust=FALSE, ...){
              if(hclust) {
                  heatmap(x@kernel, symm=TRUE, col=rev(heat.colors(n=20)), Colv=NA,labRow=NA,labCol=NA, main=list(paste('Genetic Similarity Kernel Matrix for Pathway',x@pathway@id), cex=1.4), ...)
              }else{
                  print(levelplot(x@kernel, col.regions=rev(heat.colors(n=20)),  drop.unused.levels=FALSE, scales=list(alternating=0), main=paste('Genetic Similarity Kernel Matrix for Pathway',x@pathway@id)), ...)
              }
              invisible(NULL)
          })

# S3 plot
#plot.kernel <- function(object, hclust=FALSE, ...){
#              if(hclust) {
#                  heatmap(object@kernel, symm=TRUE, col=rev(heat.colors(n=20)), Colv=NA,labRow=NA,labCol=NA, main=list(paste('Genetic Similarity Kernel Matrix for Pathway',object@pathway@id), cex=1.4), ...)
#              }else{
#                  print(levelplot(x@kernel, col.regions=rev(heat.colors(n=20)),  drop.unused.levels=FALSE, scales=list(alternating=0), main=paste('Genetic Similarity Kernel Matrix for Pathway',object@pathway@id)), ...)
#              }
#              invisible(NULL)
#          }

# kernel object constructor
setGeneric('kernel', function(object, ...) standardGeneric('kernel'))

#' \code{kernel} is a kernel object constructor and creates a kernel to be evaluated in the logistic kernel machine test. 
#'
#' @param type character, kernel type: Use \code{"lin"} for linear kernels, \code{"sia"} for size-adjusted or \code{"net"} for network-based kernels.
#' @param data Object of the class \code{GWASdata} containing the genotypes of the individuals for which a kernel will be calculated.
#' @param pathway Object of the class \code{pathway} specifying the SNP set for which a kernel will be calculated.
#'
#' @return Returns object of class kernel, including the similarity matrix of the corresponding pathway for the considered individuals (\code{kernel}).
#' @details
#' Different types of kernels can be constructed:
#' \itemize{
#'   \item \code{type='lin'} creates the linear kernel assuming additive SNP effects to be evaluated in the logistic kernel machine test. 
#'   \item \code{type='sia'} calculates the size-adjusted kernel which takes into consideration the numbers of SNPs and genes in a pathway to correct for size bias.
#'   \item \code{type='net'} calculates a network-based kernel to be evaluated in the logistic kernel machine test. Not only information on gene membership and gene/pathway size in number of SNPs is incorporated in the calculation, but also the interaction structure of genes in the pathway.
#' }
#' For more details, check the references.
#' @references
#' \itemize{
#'  \item Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin X: Powerful SNP-Set Analysis for Case-Control Genome-Wide Association Studies. Am J Hum Genet 2010, 86:929-42
#'  \item Freytag S, Bickeböller H, Amos CI, Kneib T, Schlather M: A Novel Kernel for Correcting Size Bias in the Logistic Kernel Machine Test with an Application to Rheumatoid Arthritis. Hum Hered. 2012, 74(2): 97108.
#'  \item Freytag S, Manitz J, Schlather M, Kneib T, Amos CI, Risch A, Chang-Claude J, Heinrich J, Bickeböller H: A network-based kernel machine test for the identification of risk pathways in genome-wide association studies. Hum Hered. 2013, 76(2):64-75.
#' }
#' 
#' @author Juliane Manitz, Saskia Freytag, Ngoc Thuy Ha
#' @rdname kernel-class
#' @seealso \code{\link{kernel-class}}, \code{\link{GWASdata-class}}, \code{\link{pathway-class}}
setMethod('kernel',
          definition = function(type = c('lin', 'sia', 'net'), data, pathway,
                                parallel = c('none', 'cpu', 'gpu'), ...) {
              type <- match.arg(type)
              parallel <- match.arg(parallel)
              
              if(type=='lin') k <- new('lin_kernel', data, pathway, parallel, ...)
              if(type=='sia') k <- new('sia_kernel', data, pathway, parallel, ...)
              if(type=='net') k <- new('net_kernel', data, pathway, parallel, ...)
              return(k)
})

# create linear kernel # subclass of kernel!
lin_kernel <- setClass('lin_kernel', contains = 'kernel')
# linear kernel object constructor
setGeneric('lin_kernel', function(object, ...) standardGeneric('lin_kernel'))
# @describeIn kernel
setMethod('lin_kernel', 
          definition = function(data, pathway, 
                       parallel = c('none', 'cpu', 'gpu'), ...) {
    parallel <- match.arg(parallel)
    if (inherits(data, "GWASdata")) {
        data <- data@geno }
    if (!inherits(data, "ffdf"))
        stop("not a ffdf object")
    if (is.null(attr(data, "anno")))
        stop("SNP data needs annotation as ",
             sQuote('attr(, "anno")'))

    anno <- attr(data, "anno")
    ## which SNPs are in specified pathway
    SNPset <- unique(anno$snp[which(anno$pathway == pathway@id)])
    ## subset genotype data for specified SNP set
    z <- as(data[,as.character(SNPset)],'matrix')
    if(any(is.na(z))){stop("genotype information contains missing values")}

    # K=ZZ' kernel matrix = genetic similarity
    if(parallel=='none'){
        k <- tcrossprod(z)
    }
    if(parallel=='cpu'){
        stop('sorry, not yet defined')
    }
    if(parallel=='gpu'){   ##### Suggests gputools #######
      if('gputools' %in% (.packages(all.available=TRUE))){
        require(gputools)
        z <- as.numeric(z)
        k <- gpuMatMult(z,t(z))
      }else{
        stop("Please install package 'gputools' to run matrix multiplication on GPU")
      }
    }
    k <- make_posdev(k)
    #return kernel object
    return(kernel(type='linear', kernel=k, pathway=pathway))
})

# create size-adjusted kernel
sia_kernel <- setClass('sia_kernel', contains = 'kernel')
# object constructor
setGeneric('sia_kernel', function(object, ...) standardGeneric('sia_kernel'))
#' @describeIn kernel
setMethod('sia_kernel',
          definition = function(data, pathway, 
                       parallel = c('none', 'cpu', 'gpu'), ...) {
    parallel <- match.arg(parallel)
    if (inherits(data, "GWASdata")) {
        data <- data@geno }
    if (!inherits(data, "ffdf"))
        stop("not a ffdf object")
    if (is.null(attr(data, "anno")))
        stop("SNP data needs annotation as ",
             sQuote('attr(, "anno")'))

    genemat <- function(g, data, anno){
        SNPset <- unique(anno[which(anno[,"gene"]==g),"snp"])
        z <- as(data[,as.character(SNPset)],'matrix')
        if(any(is.na(z))){stop("genotype information contains missing values")}
        z <- z[, apply(z,2,sum)/(2*nrow(z)) >= 0.001 ] #only snps maf >= 0.1%
        e.val <- eigen(cor(z), symmetric=TRUE, only.values=TRUE)$values
        nn    <- length(e.val)
        a <- matrix( rep(rowSums(z*z),nrow(z)),nrow=nrow(z))
        distances <- a -2*tcrossprod(z) + t(a)
        distances <- round(distances, digits=3)
     return( list(distances, (nn*(1-(nn-1)*var(e.val)/(nn^2))), ncol(z)) ) }
                  #matrix,         eff.length.gene,           length.gene

    genemat2 <- function(l, max.eff){
        delta <- sqrt(l[[2]]/max.eff) #delta <- sqrt(eff.length.gene/max.eff)
        roh   <- (mean(c(l[[1]])))^(-delta)*(l[[2]]/l[[3]])^(-delta)
    return(-roh*(l[[1]]/l[[3]])^(delta)) }

    anno <- attr(data, "anno")
    anno <- anno[anno[,"pathway"]==pathway@id,c("gene","snp")]#anno subset for pathway
    gene.counts <- table(anno[,"gene"]) #counts number of different SNPs per gene
    g.10 <- names(gene.counts[gene.counts >= 2]) #genes with >= 2 snps

    #[[1]]:genematrix, [[2]]:eff.length.gene, [[3]]:length.gene
    liste <- lapply(g.10, genemat, data, anno)
    get2    <- function(l){ return(l[[2]]) }
    max.eff <- max( unlist( lapply(liste, get2) ) )

    kerneltimes <- matrix( rep(0,(nrow(data))^2), nrow=nrow(data))
    kerneltimes <- Reduce('+', lapply(liste,genemat2,max.eff))
    k <- exp( sqrt(1/(length(unique(anno[,"gene"])))) * kerneltimes )
    k <- make_posdev(k)
    #return kernel object
    return(kernel(type='size-adjusted',kernel=k,pathway=pathway))
})


# network-based kernel
net_kernel <- setClass('net_kernel', contains = 'kernel')

# object constructor
setGeneric('net_kernel', function(object, ...) standardGeneric('net_kernel'))
#' @describeIn kernel
setMethod('net_kernel',
          definition = function(data, pathway, 
                       parallel = c('none', 'cpu', 'gpu'), ...) {
    if (inherits(data, "GWASdata")) {
        data <- data@geno }
    if (!inherits(data, "ffdf"))
        stop("not a ffdf object")
    if (is.null(attr(data, "anno")))
        stop("SNP data needs annotation as ",
             sQuote('attr(, "anno")'))

    anno <- attr(data, "anno")

    #genotype matrix Z, which SNPs are in specified pathway
    SNPset <- unique(anno$snp[which(anno$pathway==pathway@id)])
    #subset genotype data for specified SNP set
    Z <- as(data[,as.character(SNPset)],'matrix')
    if(any(is.na(Z))){stop("genotype information contains missing values")}

    # compute kernel
    ANA <- get.ana(anno, SNPset, pathway)
    K = Z %*% ANA %*% t(Z)

    #return kernel object
    return(kernel(type='network',kernel=K,pathway=pathway))
})

#' Apply two-step network to rewire network genes if it contains no SNPs in GWASdata (for internal use)
#'
#' @export
#' @author Juliane Manitz
#'
#' @param N adjacency matrix
#' @param remov indication which genes should be removed
#' @return An adjacency matrix containing rewired network
#'
#' @references TODO Newman?
rewire_network <- function(N, remov) {
    # early exist if no genes have to be removed
    if(length(remov)==0){ return(N) }

    # identify genes that need to be carried forward to the subnetwork
    ind_sub <- which(N[remov,] != 0)
    # extract the subnetwork
    Nsub <- N[ind_sub, ind_sub]
    # exclude self-interaction
    diag(Nsub) <- 0

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
    return(N[-remov,-remov])
}

#' Produce middle part of Network Kernel (for internal use)
#'
#' @export
#' @author Juliane Manitz, Saskia Freytag, Stefanie Freidrichs
#'
#' @param anno \code{data.frame} with annotation information
#' @param SNPset vector with SNPs to be analyzed
#' @param pathway pathway object
#' @return matrix ANA' for inner part of network kernel
get_ana <- function(anno, SNPset, pathway){

    N <- as.matrix(pathway@adj)
    N[N!=0] <- pathway@sign
    if(any(is.na(N))){stop("network information contains missing values")}

    ### remove genes that have no SNPs (not in anno):
    # genes in pathway
    net_genes <- get_genes(pathway)
    # genes in annotation -> genes with SNPs in GWASdata
    anno_sub <- anno[anno$pathway==pathway@id,]
    anno_genes <- unique(anno_sub$gene)
    # pathway genes that are not in annotation
    remov <- which(! net_genes %in% anno_genes)
    # rewire network -> separate function
    N <- rewire_network(N, remov)

    # include selfinteractions for main effects
    diag(N) <- 1
    #ensure positive definiteness by network weighting
    N <- make_posdev(N)

    #A: SNP to gene mapping
    Atab <- table(anno_sub[c('snp','gene')])
    Amat <- as(Atab, 'matrix')
    A <- Amat[SNPset,rownames(N)]    #A is colnames(Z) x rownames(N)

    #A*: size-adjustement for no of SNPs in each gene
    A.star <- A/colSums(A)

    return(A.star %*% N %*% t(A.star))
}

#' Adjust matrix to be positive definite 
#'
#' @export
#' @author Juliane Manitz, Saskia Freytag, Stefanie Freidrichs
#'
#' @param N A kernelmatrix. 
#' @return matrix N, if it was positive definite. The closes positive definite 
#' matrix to N if N was not positive definite. 
#' For more details check 
#' @references
#' \itemize{
#'  \item Freytag S, Manitz J, Schlather M, Kneib T, Amos CI, Risch A, Chang-Claude J, Heinrich J, Bickeböller H: A network-based kernel machine test for the identification of risk pathways in genome-wide association studies. Hum Hered. 2013, 76(2):64-75.
#' }
#'  
make_posdev <- function(N) {
    lambda <- min(eigen(N, only.values = TRUE, symmetric = TRUE)$values)
    # smallest eigenvalue negative = not semipositive definite
    if (lambda < 0) {
        rho <- 1/(1-lambda)
        N <- rho * N + (1-rho) * diag(dim(N)[1])
        ## now check if it is really positive definite by recursively calling
        N <- make_posdev(N)
    }
    return(N)
}
