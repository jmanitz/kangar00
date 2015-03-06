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
	if(!isSymmetric(object@kernel)){
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
              print(summary(as.vector(kern@kernel)))
              cat(paste('\nNumber of Individuals:',dim(object@kernel)[1]),'\n')
#              SNPtab <- GeneSNPsize(kernel@GWASdata)
#	      SNPno <- SNPtab[which(rownames(SNPtab)==pnet@id),]
              cat(paste('The pathway contains',dim(object@pathway@adj)[1],'genes.\n'))#,'genes and', SNPno[2], 'SNPs.\n'))
              invisible(NULL)
          })

# plot method
setGeneric('plot', function(object, ...) standardGeneric('plot'))

setMethod('plot', signature='kernel',
          definition = function(object, hclust=FALSE, ...){
		if(hclust){
		
		heatmap(object@kernel, symm=TRUE, col=rev(heat.colors(n=20)), Colv=NA,labRow=NA,labCol=NA, main=list(paste('Genetic Similarity Kernel Matrix for Pathway',object@pathway@id), cex=1.4), ...)
		}else{
		print(levelplot(object@kernel, col.regions=rev(heat.colors(n=20)),  drop.unused.levels=FALSE, scales=list(alternating=0), main=paste('Genetic Similarity Kernel Matrix for Pathway',object@pathway@id)), ...)
		}
              invisible(NULL)
          })

# create linear kernel # subclass of kernel?
kernel.lin <- function(GWASdata, pathway, parallel=c('none','cpu','gpu'), ...){
    parallel <- match.arg(parallel)
    # which SNPs are in specified pathway
    SNPset <-GWASdata@anno$rsNumber[which(GWASdata@anno$Pathway==pathway@id)]
    # subset genotype data for specified SNP set
    z <- as(GWASdata@geno[,as.character(SNPset)],'matrix') 
    # K=ZZ' ... kernel matrix = genetic similarity
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
    # return kernel object
    return(kernel(type='linear', kernel=k, pathway=pathway))
}

# create size-adjusted kernel
kernel.sia <- function(GWASdata, pathway, parallel='none', ...){
    parallel <- match.arg(parallel,c('none','cpu','gpu'))

	dat <- GWASdata@geno
	anno <- GWASdata@anno
	n <- pathway@id 

  EffectiveSNPs <- function(g.10, dat, genes){
     get.eff <- function(g,dat,genes){ 
      z <- as(dat[,as.character( genes[genes[,1]==g,2] )], "matrix")       
      z <- z[, apply(z,2,sum)/(2*dim(z)[1]) >= 0.001 ] #only snps maf >= 0.1%      
      e.val <- eigen(cor(z), symmetric=TRUE, only.values=TRUE)$values 
      return(length(e.val)*(1-(length(e.val)-1)*var(e.val)/(length(e.val)^2))) }
  return( lapply(g.10, get.eff, dat, genes) ) } 
                          
   genes <- anno[anno[,1]==paste(n,sep=""),c(2,4)] #genes and snps in pathway
   gene.counts <- table(as.character(genes[!duplicated(genes),1])) 
   #counts number of different SNPs per gene 
   g.10 <- names(gene.counts[gene.counts >= 10]) #genes with >= 10 snps 
   effSNPs <- cbind(g.10,rep(0,length(g.10)))
   effSNPs[,2] <- unlist(EffectiveSNPs(g.10, dat, genes))#, n.cores, cl)) 
   kerneltimes <- matrix( rep(0,(dim(dat)[1])^2), nrow=dim(dat)[1]) 

   g.sum <- function(g, dat, genes, effSNPs){  
      z <- as(dat[,as.character( genes[genes[,1]==g,2] )], "matrix") #rs in g       
      z <- z[, apply(z,2,sum)/(2*dim(z)[1]) >= 0.001 ] #only snps maf >= 0.1% 
      
      a <- matrix( rep(rowSums(z*z),dim(z)[1]),nrow=dim(z)[1])          
      distances <- a -2*tcrossprod(z) + t(a)
      distances <- round(distances, digits=3)
        
      length.gene <- dim(z)[2] #num. snps
      eff.length.gene<-as.numeric(effSNPs[which(as.character(effSNPs[,1])==g),2])  
      max.eff.length <- max(as.numeric(effSNPs[,2]), na.rm=T)
       
      delta <- sqrt(eff.length.gene/max.eff.length)
      roh  <- (mean(c(distances)))^(-delta)*(eff.length.gene/length.gene)^(-delta)
      return(-roh*(distances/length.gene)^(delta))  } #for one gene

  kerneltimes <- Reduce('+', lapply(g.10, g.sum, dat, genes, effSNPs)) 
  k <- exp( sqrt(1/(length(unique(genes[,1])))) * kerneltimes ) 
      
  # return kernel object
  return(kernel(type='size-adjusted',kernel=k,pathway=pathway))
} 


# network-based kernel
kernel.net <- function(GWASdata, pathway){
    source("src/get.ana.r")
	  dat <- GWASdata@geno
	  anno <- GWASdata@anno
	  n <- pathway@id
 
    #genotype matrix Z, which SNPs are in specified pathway
    SNPset <- GWASdata@anno$rsNumber[which(GWASdata@anno$Pathway==pathway@id)]
    #subset genotype data for specified SNP set 
    Z <- as(GWASdata@geno[,as.character(SNPset)],'matrix')  
    
    ANA <- get.ana(anno, SNPset, pathway)
    K = Z %*% ANA %*% t(Z)    
    #return kernel object
    return(kernel(type='network',kernel=K,pathway=pathway))
}

#function to produce middle part of net kernel:
get.ana <- function(anno, SNPset, pathway){

    N <- as.matrix(pathway@adj)
    N[N!=0] <- pathway@sign
              
    #remove genes that have no SNPs (not in anno):
    all.genes  <- colnames(N)  
    anno.genes <-unique(GWASdata@anno[which(GWASdata@anno$Pathway==pathway@id),2])
    remov <- all.genes[which(all.genes %in% anno.genes==FALSE)]      

    #rewire: may be use matrix multiplication
    if(length(remov)!=0){    
    for(g in remov){
        z <- which(colnames(N)==g)
        vec <- rbind( N[z,],seq(1:length(N[z,])) ) #column of gene to be removed
        vec <- vec[,vec[1,]!=0] #where gene has edges
        if(length(vec)>0){      #only if edges exist
          for(i in 1:(length(vec[1,])-1) ){   #i ist aktuelle edge
              for( j in (i+1):length(vec[1,])){
               if(N[vec[2,i],vec[2,j]]!=0){
		warning('Edge will be removed!')
		}
	N[vec[2,i],vec[2,j]] <- N[vec[2,i],vec[2,j]] + vec[1,i]*vec[1,j]
               N[vec[2,j],vec[2,i]] <- N[vec[2,j],vec[2,i]] + vec[1,i]*vec[1,j]}  
          }    #additive
        }} 
     N <- N[-z,-z] }
     
    #check, if both directions exist and have same relation (value 2/-2):
    N[N>1]    <-  1 
    N[N<(-1)] <- -1    
    # include main effects
    diag(N) <- 1  
    
    #ensure positive definiteness by network weighting
    lambda <- min(eigen(N, only.values=TRUE, symmetric=TRUE)$values)
    if(lambda<0){ # smallest eigenvalue negative = not semipositive definite
        rho <- 1/(1-lambda)
        N <- rho*N + (1-rho)*diag(dim(N)[1])
    }    
          
    #A: SNP to gene mapping      
    A.table <- t(table(anno))
    A <- A.table[SNPset,rownames(N)]    #= A.table[colnames(Z),rownames(N)]

    #A*: size-adjustement for no of SNPs in each gene
    A.star <- A/colSums(A)
    
return(A.star %*% N %*% t(A.star))
}
