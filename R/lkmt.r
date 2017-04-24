#################################################
#
# logistic kernel maschine test functions
#
#################################################


#' An S4 class to represent the variance component test.
#'
#' @slot formula A formula stating the regression nullmodel that will be used in
#' the variance component test.
#' @slot kernel An object of class \code{\link{kernel}} representing the similarity
#' matrix of the individuals based on which the pathways influence is evaluated.
#' @slot GWASdata An object of class \code{\link{GWASdata}} including the data
#' on which the test is conducted.
#' @slot statistic A \code{vector} giving the value of the variance component
#' test statistic.
#' @slot df A \code{vector} containing the number of degrees of freedom.
#' @slot p.value A \code{vector} giving the p-value calculated for the pathway
#' in the variance component test.
#'
#' For details on the variance component test see the references.
#' @references
#' \itemize{
#' \item  Liu D, Lin X, Ghosh D: Semiparametric regression of multidimensional genetic pathway data: least-squares kernel machines and linear mixed models. Biometrics 2007, 63(4):1079-88.
#'  \item Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin X: Powerful SNP-Set Analysis for Case-Control Genome-Wide Association Studies. Am J Hum Genet 2010, 86:929-42
#' }
#'
#' @author Juliane Manitz, Stefanie Friedrichs
#' @rdname lkmt-class
#' @export lkmt
#' @import methods
lkmt <- setClass('lkmt',
                 slots=c(formula='formula', kernel='kernel', GWASdata='GWASdata',
                         statistic='vector',df='vector',p.value='vector'))

setValidity('lkmt', function(object){ 
 	msg  <- NULL
	valid <- TRUE
  if(!(class(object@GWASdata)=="GWASdata")){
	  valid=FALSE
	  msg <- c(msg, "no GWASdata object!")
	}
 	if(!is.matrix(object@kernel)){
	  valid=FALSE
	  msg <- c(msg, "kernel matrix is missing!")
	}
  if(!is.numeric(object@statistic)){
	  valid=FALSE
	  msg <- c(msg, "test statistic value must be numeric!")
	}
  if(!is.numeric(object@df)){
	  valid=FALSE
	  msg <- c(msg, "degrees of freedom must be given!")
	}
  if(!is.numeric(object@p.value)){
	  valid=FALSE
	  msg <- c(msg, "p-value needs to be numeric!")
	}
})

#' \code{show} Shows basic information on \code{lkmt} object
#'
#' @param object An object of class \code{lkmt}.
#' @return \code{show} Basic information on \code{lkmt} object.
#' @examples
#' # show method
#' data(lkmt.net.kernel.hsa04020)
#' lkmt.net.kernel.hsa04020
## @author Juliane Manitz
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
## @param object An object of class \code{\link{lkmt}}.
#' @param ... Further arguments can be added to the function
#' @return \code{summary} Summarized information on \code{lkmt} object.
#' @examples
#' # summary method
#' summary(lkmt.net.kernel.hsa04020)
## @author Juliane Manitz
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


#' A function to calculate the p-values for kernel-matrices. 
#'
#' @param formula The formula to be used for the regression nullmodel.  
#' @param kernel An object of class \code{kernel} including the pathway 
#' representing kernel-matrix based on which the test statistic will be calculated.
#' @param GWASdata A \code{GWASdata} object stating the data used in analysis. 
#' @param method A \code{character} specifying which method will be used for 
#' p-value calculation. Available are \code{'satt'} for the Satterthwaite 
#' approximation and \code{'davies'} for Davies' algorithm. For more details 
#' see the references.
#' @param ... Further arguments can be given to the function.
#' @return 
#' An \code{lkmt} object including the following test results
#' \itemize{
#' \item The formula of the regression nullmodel used in the variance 
#' component test.
#' \item An object of class \code{\link{kernel}} including the similarity
#' matrix of the individuals based on which the pathways influence is evaluated.
#' \item An object of class \code{\link{GWASdata}} stating the data on
#' which the test was conducted.
#' \item statistic A \code{vector} giving the value of the variance component
#' test statistic.
#' \item df A \code{vector} giving the number of degrees of freedom.
#' \item p.value A \code{vector} giving the p-value calculated for the pathway
#' in the variance component test.
#' }
#' @references
#' For details on the variance component test 
#' \itemize{
#'  \item Wu MC, Kraft P, Epstein MP, Taylor DM, Chanock SJ, Hunter DJ, Lin X: Powerful SNP-Set Analysis for Case-Control Genome-Wide Association Studies. Am J Hum Genet 2010, 86:929-42
#' \item  Liu D, Lin X, Ghosh D: Semiparametric regression of multidimensional genetic pathway data: least-squares kernel machines and linear mixed models. Biometrics 2007, 63(4):1079-88.
#' }
#' @examples
#' data(net.kernel.hsa04020)
#' data(gwas)
#' lkmt_test(pheno ~ sex + age, net.kernel.hsa04020, gwas, method='satt')
#' @author Stefanie Friedrichs, Juliane Manitz
#' @export
#' @rdname lkmt_test
lkmt_test <- function(formula, kernel, GWASdata, method=c('satt','davies'), ...){

    if(length(GWASdata@pheno) == 0) stop("Please specify phenotypes.")

    method <- match.arg(method)
    nullmodel <- stats::glm(formula, data=GWASdata@pheno, family=stats::binomial, x=TRUE)
    if(method == 'satt'){
        model <- score_test(kernel@kernel, nullmodel)
    }
    if(method == 'davies'){
        model <- davies_test(kernel@kernel, nullmodel)
    }
    ret <- new('lkmt', formula=formula, kernel=kernel, GWASdata=GWASdata,
               statistic=model$statistic, df=model$parameter, p.value=model$p.value)
    return(ret)
}

setGeneric('score_test', function(x1, x2, ...) standardGeneric('score_test'))
#' Calculates the p-value for a kernelmatrix using Satterthwaite approximation
#'
#' For parameter \code{"satt"} a pathways influence on the probability of 
#' beeing a case is evaluated in the logistic kernel machine test and p-values 
#' are determined using a Sattherthwaite Approximation as described by Dan Schaid.
#'
#' @param x1 A \code{\link{matrix}} which is the
#' similarity matrix calculated for the pathway to be tested.
#' @param x2 An \code{lm} or \code{glm} object of the nullmodel with fixed 
#' effects covariates included, but no genetic random effects.
## @return A \code{list} including the test results of the pathway.
## @author Stefanie Friedrichs, Saskia Freytag, Ngoc-Thuy Ha
#' @references
#' For details on the p-value calculation see
#' \itemize{
#' \item Schaid DJ: Genomic Similarity and Kernel Methods I: Advancements by 
#' Building on Mathematical and Statistical Foundations. Hum Hered 2010, 70:109-31
#' }
#' @rdname lkmt_test
setMethod('score_test', signature(x1 = 'matrix'), 
          definition = function(x1, x2){   
        if(sum(is(x2) %in% c("glm","lm")) == 0 ){
                stop("nullmodel should be a glm- or lm-object!")
        }
        if(is.null(x2$x)){
                stop("The glm-object should have a design-matrix x!")
        }
        nas <- x2$na.action
        Y   <- x2$y
        X   <- x2$x
        mui <- x2$fitted.values
        W <- (mui*(1-mui)) * diag(length(mui))
        #P <- W-W%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W
        WX <- W%*%X
        P  <- W-WX%*%solve(t(X)%*% WX, t(X)%*%W)
        if(!is.null(nas)){ x1 <- x1[-nas,-nas] }
        
        TT <- 1/2*t(Y-mui)%*%x1%*%(Y-mui)
        PK <- P%*%x1
        A  <- W%*% PK %*%P%*%W
        VarT <- 1/2*sum(diag((PK%*%PK)))+
                1/4*sum(diag(A)^2*mui*(1-mui)*(1+6*mui^2-6*mui))
        ExpT <- 1/2*sum(diag((PK)))

        a <- VarT/(2*ExpT)
        d <- (2*ExpT^2)/VarT
        p.value <- stats::pchisq(TT/a, d, lower.tail=FALSE)

        mod <- list(p.value=p.value,          
               statistic=as.numeric(TT/a),  
               parameter=d,                 
               method="Dan Schaid Score Test") 
        return(mod) 
})

setGeneric('davies_test', function(x1, x2, ...) standardGeneric('davies_test'))
#' Calculates the p-value for a kernel matrix using davies method
#'
#' For parameter \code{"davies"} a pathways influence on the probability
#' of beeing a case is evaluated using the p-value calculation method described 
#' by Davies. Here the function \code{davies()} from package \code{CompQuadForm}
#' is used.
#'
## @param x1 A \code{matrix} which is the
## similarity matrix calculated for the pathway to be tested.
## @param x2 A \code{glm} object of the nullmodel with fixed effects
## covariates included, but no genetic random effects.
## @return A \code{list} including the test results of the pathway.
## @author Stefanie Friedrichs
#' @references
#' \itemize{
#' \item Davies R: Algorithm as 155: the distribution of a linear combination of
#'      chi-2 random variables. J R Stat Soc Ser C 1980, 29:323-333.
#' }
#' @import CompQuadForm
#' @rdname lkmt_test
setMethod('davies_test', signature(x1 = 'matrix'), 
          definition = function(x1, x2){
          
        if(sum(is(x2) %in% c("glm","lm")) == 0 ){
                stop("nullmodel should be a glm- or lm-object!")
        }
        if(is.null(x2$x)){
                stop("The glm-object should have a design-matrix x!")
        }

        nas <- x2$na.action
        if(!is.null(nas)){
           x1 <- x1[-nas,-nas]
        }

        Y   <- x2$y
        X   <- x2$x                  
        mui <- x2$fitted.values

        if( min(eigen(x1)$values)<0 ){ x1 <- make_psd(x1) }

        Q <- 1/2*t(Y-mui)%*%x1%*%(Y-mui)
        P <- diag(nrow(X))-X%*%solve(t(X)%*%X,t(X))#solve(t(X)%*%X)%*%t(X)

        W <- (mui*(1-mui))*diag(length(mui))
        v <- eigen(W)$vectors  
        Wsqrt <- v %*% diag(sqrt(eigen(W)$values)) %*% t(v)
        WP    <- Wsqrt%*%P

        lambda <- eigen(0.5*WP%*%x1%*%t(WP))$values
        #delta <- rep(0, length(lambda))
        #sigma <- 0
        pval <- davies(Q,lambda, rep(1, length(lambda)), rep(0, length(lambda)),
                       0, 10000, 0.0001)$Qq  #needs CompQuadForm

        #p-value of test, test staistic value, degree of freedom, name
        mod <- list(p.value=pval, statistic=Q, parameter=1,
                   method="Davies Score Test")
        names(mod$statistic) <- "Chi Squared Test Statictic Value"
        names(mod$parameter) <- "Degree of Freedom"

   return(mod)
})
