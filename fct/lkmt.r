#################################################
#
# loistic kernel maschine test functions
#
#################################################

# object constructor 
#!! define lkmt object !!
lkmt <- setClass('lkmt',
                  slots=c(formula='formula', kernel='kernel', GWASdata='GWASdata',statistic='vector',df='vector',p.value='vector'))

# type ... lin, adj, or net kernel
# kernel .. kernel matrix of dimension equal to individuals
# GWASdata .. GWASdata object
# pathway .. pathway information

setValidity('lkmt', function(object){  # to be defined !!
	msg  <- NULL
	valid <- TRUE
	print('create validity function!')
	if(valid) TRUE else msg
})

lkmt <- function(formula, kernel, GWASdata, ...){
    nullmodel<-glm(formula, data=GWASdata@pheno, family=binomial, x=TRUE)
    model <- score.test(kernel@kernel, nullmodel, pd.check=FALSE)[[1]]
    # return object
    ret <- new('lkmt', formula=formula, kernel=kernel, GWASdata=GWASdata, statistic=model$statistic, df=model$parameter, p.value=model$p.value)
    return(ret)
}

# show method
setMethod('show', signature='lkmt',
          definition = function(object){
              cat('An object of class ', class(object), ':\n', sep='')
              cat('GWAS data:', object@GWASdata@desc,'\n')
              cat('Kernel:', object@kernel@type,'\n')
              cat('Pathway:', object@kernel@pathway@id,'\n\n')
	      cat('The score test results in the p-value: ',object@p.value,'\n\n')
              invisible(NULL)
          })

# summary method
setGeneric('summary', function(object, ...) standardGeneric('summary'))

setMethod('summary', signature='lkmt',
          definition = function(object){
              cat('An object of class ', class(object), ':\n\n', sep='')
              cat('GWAS data:', object@GWASdata@desc,'\n')
              cat('Kernel:', object@kernel@type,'\n')
              cat('Pathway:', object@kernel@pathway@id,'\n\n')
	      cat('\tDan Schaid Score Test\n\n')
	      cat('Chi Squared Test Statictic Value =', object@statistic,'\nDegrees of Freedom =',object@df,'\np-value =', object@p.value,'\n\n')
#	      cat('The score test results in the p-value: ',object@pval,'\n\n')
              invisible(NULL)
          })

############################################################
###    Score Test Using Sattherthwaite Approximation     ###             
###                ( Dan Schaid )                        ###             
############################################################

#require(Matrix)

score.test <- function(
                       kernels,        # list of kernels-matrices, or a kernel-matrix
                       nullmodel,      # a glm-object
                       pd.check=TRUE,  # boolean, whether to check for positive definiteness
                       kernel.out=FALSE# boolean, whether the kernel need to be omitted
                       ){

        ## check whether kernel is matrix or list
        if(is.matrix(kernels)){
                kernels <- list(pathway1=kernels)
        }else if(!is.list(kernels)){
                stop("kernels should be either a kernel-matrix or a list of kernel matrices!")
        }
        
        ## check whether nullmodel is glm-object
        if(sum(is(nullmodel) %in% c("glm","lm")) == 0 ){
                stop("nullmodel should be a glm- or lm-object!")
        }
        
        ## check whether nullmodel has x
        if(is.null(nullmodel$x)){
                stop("The glm-object should have a design-matrix x!")
        }
                
        nas <- nullmodel$na.action
        
                        
        pathwaynames <- names(kernels)
        Y <- nullmodel$y
        X <- nullmodel$x
        mui <-  nullmodel$fitted.values
        W <- Diagonal(length(mui),mui*(1-mui))
        #P <- W-W%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W 
        WX <- W%*%X
        P <- W-WX%*%solve(t(X)%*% WX, t(X)%*%W) 
        
        all.mod <- list()
        for(k in pathwaynames){

                K <- kernels[[k]]
                K <- as.matrix(K)
                if(!is.null(nas)){
                        K <- K[-nas,-nas]
                }
                
                ### check, whether K is positive definit:
                #try1 <- try(chol(K))
                #if(is(try1) == "try-error"){
                #        if(pd.check==T){
                #        eigval <- eigen(K)$values
                #        if(min(eigval,na.rm=T) < 0){
                #                warning(paste("The Kernel of pathway", k, "may be not semi-positive definit! Minimal eigenvalue: ", min(eigval,na.rm=T)))
                #        }
                #    }
                #} 
                           
                      
                TT <- 1/2*t(Y-mui)%*%K%*%(Y-mui)
                PK <- P%*%K
                A <- W%*% PK %*%P%*%W  
                VarT <- 1/2*sum(diag((PK%*%PK)))+
                        1/4*sum(diag(A)^2*mui*(1-mui)*(1+6*mui^2-6*mui))
                ExpT <- 1/2*sum(diag((PK)))


                a <- VarT/(2*ExpT)
                d <- (2*ExpT^2)/VarT
                p.value <- pchisq(TT/a, d, lower.tail=FALSE)
        
                mod <- list(p.value=p.value,    # p-value of test
                        statistic=TT/a,         # the test staistic value
                        parameter=d,             # the degree of freedom
                        method="Dan Schaid Score Test" # name of test
                        )
                names(mod$statistic) <- "Chi Squared Test Statictic Value"
                names(mod$parameter) <- "Degree of Freedom"
                if(kernel.out==T){
                        mod$kernel <- K         # if needed, the kernel wil be also omitted
                }
                class(mod) <- "htest"
        
                all.mod[[as.character(k)]] <- mod
        } # end of for pathwaynames
        return(all.mod)
}

