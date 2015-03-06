#################################################
#
# GWASdata object functions
#
#################################################

# object constructor
GWASdata <- setClass('GWASdata', slots=c(pheno='data.frame', geno='databel', anno='data.frame', desc='character'))

# pheno .. data.frame specifying ids, phenotypes and covariates e.g. ID, pheno, sex, pack.years
# geno .. genotype information in databel format
# anno .. data.frame; annotation file mapping SNP -> genes -> pathways
# desc ... character giving GWAS description and information

# validy checks
setValidity('GWASdata', function(object){  # to be defined !!
	msg  <- NULL
	valid <- TRUE
# validy=check.GWASdata(dat=geno,anno=anno,ids=rownames(pheno)) # to be defined !!
#check.GWASdata <- function(dat,anno,ids){ 
# !! this function needs improvement !!

#        test1 <- try(databel(dat))
#        if(is(test1) == "try-error"){
#                stop("Did not find databel data!")
#        }
        if(!is.data.frame(object@anno)){
		valid <- FALSE
                msg <- c(msg, "anno need to be a data frame!")
        }
#	ids <- rownames(object@pheno)
#        if(!("vector" %in% is(ids))){
#		valid <- FALSE
#                msg <- c(msg, "pathway.ids need to be a vector of pathway ID's!")
#        }
#        test2 <- ids %in% anno[,1]
#        if(sum(test2) != length(test2)){
#                stop("There are pathway ID's not in the Annotation file!")
#        }
#        disconnect(test1)
	if(valid) TRUE else msg
})

# show method
setMethod('show', signature='GWASdata',
          definition = function(object){
              cat('An object of class ', class(object), ' from ',object@desc,'\n\n',sep='')
              cat('Phenotypes for ', dim(object@pheno)[1], ' individuals: \n',sep='')
              print(head(object@pheno))
              invisible(NULL)
          })

# summary method
setGeneric('summary', function(object, ...) standardGeneric('summary'))

setMethod('summary', signature='GWASdata',
          definition = function(object){
              cat('An object of class ', class(object), ' from ',object@desc,'\n\n',sep='')
              cat('Phenotypes for ', dim(object@pheno)[1], ' individuals: \n\n',sep='')
              print(summary(object@pheno)) # summary of phenotype and covariate data
              cat('\nGenotypes: \n\n')
	      cat('Total number of genes and SNPs per pathway:\n')
              print(GeneSNPsize(object)) # overview of pathway - gene - snp mapping
 	      cat('\n')
	      print(object@geno)
              invisible(NULL)
          })

## GeneSNPsize() to be called in summary(GWASdata)
setGeneric('GeneSNPsize', function(object, ...) standardGeneric('GeneSNPsize'))

setMethod('GeneSNPsize', signature='GWASdata',
# Create list of pathway names + snp size and genes size
# counts number of Genes/SPNs in each pathway 
          definition <- function(object){
              anno <- object@anno    #Annotation file
              anno <- anno[!duplicated(anno[,1:4]),]
              nrsnps <- table(anno[,1])
              anno <- anno[!duplicated(anno[,1:2]),]
              nrgenes <- table(anno[,1])
              tab <- cbind(nrgenes,nrsnps)
	      colnames(tab) <- c('genes','SNPs')
              return(tab) 
          })

