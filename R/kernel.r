#################################################
#
# kernel object functions
#
#################################################

# object constructor
# !! integrate update function from make.BioPax.r !!
kernel <- setClass('kernel',
                   slots=c(type='character', kernel='matrix', pathway='pathway'))

# type ... lin, adj, or net kernel
# kernel .. kernel matrix of dimension equal to individuals
# GWASdata .. GWASdata object
# pathway .. pathway information

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
    #          msg <- c(msg, "kernel matrix has to be dimension equal to individuals in GWAS")
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
setMethod('show', signature='kernel',
          definition = function(object){
              cat('An object of class ', class(object), ' of type ', object@type,' for pathway ', object@pathway@id, '\n\n',sep='')
              print(object@kernel)
              invisible(NULL)
          })

# summary method
setGeneric('summary', function(object, ...) standardGeneric('summary'))

setMethod('summary', signature='kernel',
          definition = function(object){
              cat('An object of class ', class(object), ' of type ', object@type,' for pathway ', object@pathway@id, ' with values: \n\n',sep='')
              cat(paste('\nNumber of Individuals:',dim(object@kernel)[1]),'\n')
              #              SNPtab <- GeneSNPsize(kernel@GWASdata)
              #	      SNPno <- SNPtab[which(rownames(SNPtab)==pnet@id),]
              cat(paste('The pathway contains',dim(object@pathway@adj)[1],'genes.\n'))#,'genes and', SNPno[2], 'SNPs.\n'))
              invisible(NULL)
          })

# plot method
if (!isGeneric("plot")) setGeneric('plot')
#setGeneric('plot', function(object, ...) standardGeneric('plot'))
#
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


# create linear kernel # subclass of kernel?
kernel.lin <- function(data, pathway, parallel=c('none','cpu','gpu'), ...) {
    parallel <- match.arg(parallel)
    if (inherits(data, "GWASdata")) {
        data <- data@geno }
    if (!inherits(data, "databel"))
        stop("not a databel object")
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
    if(parallel=='gpu'){
        z <- magma(z, gpu=TRUE)
        k <- tcrossprod(z)
    }
    
    k <- make_posdev(k)
    #return kernel object
    return(kernel(type='linear', kernel=k, pathway=pathway))
}



# create size-adjusted kernel
kernel.sia <- function(data, pathway, parallel='none', ...){
    parallel <- match.arg(parallel,c('none','cpu','gpu'))
    if (inherits(data, "GWASdata")) {
        data <- data@geno }
    if (!inherits(data, "databel"))
        stop("not a databel object")
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
}


# network-based kernel
kernel.net <- function(data, pathway) {   
    if (inherits(data, "GWASdata")) {
        data <- data@geno }
    if (!inherits(data, "databel"))
        stop("not a databel object")
    if (is.null(attr(data, "anno")))
        stop("SNP data needs annotation as ",
             sQuote('attr(, "anno")'))
    
    anno <- attr(data, "anno")
    
    #genotype matrix Z, which SNPs are in specified pathway
    SNPset <- unique(anno$snp[which(anno$pathway==pathway@id)])
    #subset genotype data for specified SNP set
    Z <- as(data[,as.character(SNPset)],'matrix')
    if(any(is.na(Z))){stop("genotype information contains missing values")}
    
    ANA <- get.ana(anno, SNPset, pathway)
    K = Z %*% ANA %*% t(Z)
    #return kernel object
    return(kernel(type='network',kernel=K,pathway=pathway))
}

#function to produce middle part of net kernel:
get.ana <- function(anno, SNPset, pathway){
    
    N <- as.matrix(pathway@adj)
    N[N!=0] <- pathway@sign
    if(any(is.na(N))){stop("network information contains missing values")}
    #remove genes that have no SNPs (not in anno):
    all.genes  <- colnames(N)
    anno.genes <- unique(anno$gene[which(anno$pathway==pathway@id)])
    remov <- all.genes[which(all.genes %in% anno.genes==FALSE)]
    
    #rewire: may be use matrix multiplication
    if(length(remov)!=0){
    for(g in remov){   
      z <- which(colnames(N)==g)
      vec <- rbind( N[z,],seq(1:length(N[z,])) )    
      vec <- data.frame( vec[,vec[1,]!=0]) 
        #if something must be rewired       
        if( length(vec[1,])>1 ){      
          for(i in 1:(length(vec[1,])-1) ){  
              if((i+1)<=length(vec[1,])){
              for( j in (i+1):length(vec[1,]) ){ #i ist aktuelle edge
               if(N[vec[2,i],vec[2,j]]!=0){print("Edge will be removed!")}
               N[vec[2,i],vec[2,j]] <- N[vec[2,i],vec[2,j]] + vec[1,i]*vec[1,j]
               N[vec[2,j],vec[2,i]] <- N[vec[2,j],vec[2,i]] + vec[1,i]*vec[1,j]}}       
         } }                                   
     N <- N[-z,-z]
     N[N>1]    <-  1
     N[N<(-1)] <- -1  
     } }
     
    #check, if both directions exist and have same relation (value 2/-2):
    N[N>1]    <-  1
    N[N<(-1)] <- -1 
    #include main effects
    diag(N) <- 1   
    #ensure positive definiteness by network weighting
    N <- make_posdev(N)
    
    #A: SNP to gene mapping
    A.table <- t(table(anno[which(anno$pathway==pathway@id),c("gene","snp")]))
    A.table <- as(A.table, 'matrix') 
    A <- A.table[SNPset,rownames(N)]    #A is colnames(Z) x rownames(N)
     
    #A*: size-adjustement for no of SNPs in each gene
    #A.star <- A
    #for(i in 1:length(A[1,])){A.star[,i] <- A[,i]/sum(A[,i])}
    A.star <- A/colSums(A)
   
    return(A.star %*% N %*% t(A.star))
}

make_posdev <- function(N) {
    lambda <- min(eigen(N, only.values = TRUE, symmetric = TRUE)$values)
    # smallest eigenvalue negative = not semipositive definite
    if (lambda < 0) {
        rho <- 1/(1-lambda)
        N <- rho * N + (1-rho) * diag(dim(N)[1])
    }
    return(N)
}
