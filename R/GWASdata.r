#################################################
#
# GWASdata object functions
#
#################################################

# object constructor
GWASdata <- setClass('GWASdata', slots=c(pheno='data.frame', geno='databel', desc='character'))

# pheno .. data.frame specifying ids, phenotypes and covariates e.g. ID, pheno, sex, pack.years
# geno .. genotype information in databel format
# anno .. annotation file mapping SNP -> genes -> pathways
#         data.frame with variables:
#             pathway ... pathway id
#             gene    ... gene name (which version?)
#             chr     ... chromosome
#             snp     ... rsNumber
#             position... SNP position
# desc ... character giving GWAS description and information

#GWASdata@geno
#geno <- setClass('geno', slots='databel')
#attr(geno, 'anno') <- anno # data.frame

# validy checks
setValidity('GWASdata', function(object){  # !! this function needs improvement !!
	msg  <- NULL
	valid <- TRUE

        # check annotation file
	if(is.null(attr(object@geno,'anno'))){ #does an additional attribute with annotation exists for databel object?
		valid <- FALSE
                msg <- c(msg, "databel object geno needs an additional attribute: data.frame 'anno'")
	}
	if(!is.data.frame(attr(object@geno,'anno'))){
		valid <- FALSE
                msg <- c(msg, "geno attribute anno need to be a data frame")
	}
        anno.names <- c('pathway','gene','chr','snp','position')
        # check anno data frame variable names
        if(!all(anno.names %in% colnames(anno))){
		valid <- FALSE
                msg <- c(msg, paste("anno variable names to not match:",anno.names))
                anno <- anno[,match(colnames(anno), anno.names)]  
        }
# 	 # check whether GWASdata@geno has missings or not !todo!
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
              anno <- attr(object@geno, 'anno')    #Annotation file
#              anno <- anno[!duplicated(anno[,1:4]),]
#              nrsnps <- table(anno[,1])
              nrsnps <- table(unique(anno[,c('pathway','gene','snp')])$pathway)
#              anno <- anno[!duplicated(anno[,1:2]),]
#              nrgenes <- table(anno[,1])
              nrgenes <- table(unique(anno[,c('pathway','gene')])$pathway)
              tab <- cbind(nrgenes,nrsnps)
	      colnames(tab) <- c('genes','SNPs')
              return(tab) 
          })

