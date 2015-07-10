#################################################
#
# logistic kernel maschine test functions
#
#################################################


#' An S4 class to represent the variance component test. 
#' 
#' @rdname lkmt-class
#' @slot formula A formular stating the regression nullmodel that will be used in
#' the variance component test. 
#' @slot kernel An object of class \code{\link{kernel}} representing the similarity
#' matrix of the individuals based on which the pathways influence is evaluated.
#' @slot GWASdata An object of class \code{\link{GWASdata}} stating the data on 
#' which the test is conducted. 
#' @slot statistic A \code{vector} giving the value of the variance component 
#' ttest statistic.
#' @slot statistic A \code{vector} giving the number of degrees of freedom. 
#' @slot statistic A \code{vector} giving the p-value calculated for the pathway 
#' in the variance component test.
#'
#' @author Juliane Manitz, Stefanie Friedrichs
#' @export lkmt
lkmt <- setClass('lkmt',
                 slots=c(formula='formula', kernel='kernel', GWASdata='GWASdata',
                         statistic='vector',df='vector',p.value='vector'))

setValidity('lkmt', function(object){  # to be defined !!
	msg  <- NULL
	valid <- TRUE
	print('create validity function!')
	if(valid) TRUE else msg
})


#' Function to start the logistic kernel machine test. 
#' 
#' @param formula A formular for the regression nullmodel that will be used in
#' the variance component test. 
#' @param kernel An object of class \code{\link{kernel}} representing the similarity
#' matrix of the individuals based on which the pathways influence is evaluated.
#' @param GWASdata An object of class \code{\link{GWASdata}} representing the data on 
#' which the test is conducted. 
#' @param ... additional arguments can be added.
#' @return An \code{\link{lkmt}} object with the test results. 
#' @author Juliane Manitz, Stefanie Friedrichs
#'
#' For details on the variance component test see the references.
#' @references
#' \itemize{
#'  \item Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin X: Powerful SNP-Set Analysis for Case-Control Genome-Wide Association Studies. Am J Hum Genet 2010, 86:929-42
#' }
lkmt <- function(formula, kernel, GWASdata, ...){
    nullmodel <- glm(formula, data=GWASdata@pheno, family=binomial, x=TRUE)
    model     <- score_test(kernel@kernel, nullmodel, pd.check=FALSE)[[1]]
    ret       <- new('lkmt', formula=formula, kernel=kernel, GWASdata=GWASdata, 
                     statistic=model$statistic, df=model$parameter, 
                     p.value=model$p.value)
    return(ret)
}

#' \code{show} Shows basic information on \code{lkmt} object
#' 
#' @param object An object of class \code{\link{lkmt}}.
#' @return \code{show} Basic information on \code{lkmt} object.
#' #@author Juliane Manitz
#' @export
#' @rdname lkmt-class
#' @aliases show,GWASdata,ANY-method
setMethod('show', signature='lkmt',
          definition = function(object){
              cat('An object of class ', class(object), ':\n', sep='')
              cat('GWAS data:', object@GWASdata@desc,'\n')
              cat('Kernel:', object@kernel@type,'\n')
              cat('Pathway:', object@kernel@pathway@id,'\n\n')
	      cat('The score test results in the p-value: ',object@p.value,'\n\n')
              invisible(NULL)
          })

## summary
setGeneric('summary', function(object, ...) standardGeneric('summary'))

#' \code{summary} Summarizes information on \code{lkmt} object
#' 
#' @param object An object of class \code{\link{lkmt}}.
#' @return \code{summary} Summarized information on \code{lkmt} object.
#' #@author Juliane Manitz
#' @export
#' @rdname lkmt-class
#' @aliases summary,GWASdata,ANY-method
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


#' Calculates the p-value for a kernelmatrix  
#'
#' This function evaluates a pathways influence on an individuals probability
#' of beeing a case using the logistic kernel machine test. P-values are 
#' determined using a Sattherthwaite Approximation as described by Dan Schaid. 
#'
#' @param kernels an object of \code{\link{kernel}} class. Includes the 
#' similarity matrix calculated for the pathway to be tested. 
#' @param nullmodel A \code{glm} object of the nullmodel with fixed effects 
#' covariates included, but no genetic random effects. 
#' @param  pd.check boolean, whether to check for positive definiteness.
#' @return A \code{list} including the test results of the pathway.   
#' @author Stefanie Friedrichs, Saskia Freytag, Ngoc-Thuy Ha
#'
#' For details on the p-value approximation see
#' \itemize{
#' \item Schaid DJ: Genomic Similarity and Kernel Methods I: Advancements by Building on Mathematical and Statistical Foundations. Hum Hered 2010, 70:109-31
#' }

score_test <- function(kernels, nullmodel, pd.check=TRUE){ 
                                                          
        if(is.matrix(kernels)){
                kernels <- list(pathway1=kernels)
        }else if(!is.list(kernels)){
                stop("kernels should be a kernel-matrix!")
        }
        if(sum(is(nullmodel) %in% c("glm","lm")) == 0 ){
                stop("nullmodel should be a glm- or lm-object!")
        }       
        if(is.null(nullmodel$x)){
                stop("The glm-object should have a design-matrix x!")
        }
        nas <- nullmodel$na.action
        pathwaynames <- names(kernels)
        Y <- nullmodel$y
        X <- nullmodel$x
        mui <-  nullmodel$fitted.values
        #W <- Diagonal(length(mui),mui*(1-mui))  , library(Matrix)
        W <- (mui*(1-mui)) * diag(length(mui))
        #P <- W-W%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W 
        WX <- W%*%X
        P  <- W-WX%*%solve(t(X)%*% WX, t(X)%*%W) 
        all_mod <- list()
        for(k in pathwaynames){
                K <- kernels[[k]]
                K <- as.matrix(K)
                if(!is.null(nas)){
                   K <- K[-nas,-nas]
                }                                          
                TT <- 1/2*t(Y-mui)%*%K%*%(Y-mui)
                PK <- P%*%K
                A  <- W%*% PK %*%P%*%W  
                VarT <- 1/2*sum(diag((PK%*%PK)))+
                        1/4*sum(diag(A)^2*mui*(1-mui)*(1+6*mui^2-6*mui))
                ExpT <- 1/2*sum(diag((PK)))

                a <- VarT/(2*ExpT)
                d <- (2*ExpT^2)/VarT
                p.value <- pchisq(TT/a, d, lower.tail=FALSE)
        
                mod <- list(p.value=p.value,     # p-value of test
                        statistic=TT/a,          # the test staistic value
                        parameter=d,             # the degree of freedom
                        method="Dan Schaid Score Test") # name of test
                names(mod$statistic) <- "Chi Squared Test Statictic Value"
                names(mod$parameter) <- "Degree of Freedom"        
                all_mod[[as.character(k)]] <- mod
        } # end of for pathwaynames
        return(all_mod)
}

