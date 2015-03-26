#################################################
#
# GWASdata object functions
#
#################################################

# object constructor
GWASdata <- setClass('GWASdata', slots=c(pheno='data.frame', geno='databel', desc='character'))

# pheno .. data.frame specifying ids, phenotypes and covariates e.g. ID, pheno, sex, pack.years
#          Ids have to be in first column!
# geno .. genotype information in databel format
# anno .. annotation file mapping SNP -> genes -> pathways
#         data.frame with variables:
#             pathway ... pathway id
#             gene    ... gene name (which version?)
#             chr     ... chromosome
#             snp     ... rsNumber
#             position... SNP position
# desc ... character giving GWAS description and information

# validy checks
setValidity('GWASdata', function(object){
    msg   <- NULL
    valid <- TRUE

    ## check annotation file
    if( is.null( attr(object@geno, "anno") ) ){ #does attribute for databel object exists?
        valid <- FALSE
        msg   <- c(msg, "databel object geno needs an additional attribute: data.frame 'anno'")
    }
    if(!is.data.frame(attr(object@geno,"anno"))){ #is attribute dataframe?
  	valid <- FALSE
        msg   <- c(msg, "geno attribute anno needs to be a data frame")
    }
    ## check anno data frame variable names
    anno.names <- c('pathway','gene','chr','snp','position')
    if(!all(anno.names %in% colnames(attr(object@geno,'anno')))){
        valid <- FALSE
        msg   <- c(msg, paste("anno variable names do not match:",anno.names))
    }

    ## check whether GWASdata@geno has missings
    # <FIXME> Currently not working
    # if(sum(is.na(object@geno))>0){
    #     valid <- FALSE
    #     msg   <- c(msg, "genotypes include missing values and need to be imputed!")
    # }

    ## phenotypes for more individuals than have genotypes
    if(length(unique(object@pheno[,1]))>length(rownames(object@geno))){
        valid <- FALSE
        msg <- c(msg, "phenotypes exist for more individuals than have genotypes!")
    }
    ## check order of individuals in genotypes and phenotypes
    if(!all.equal(as.character(object@pheno[,1]),rownames(object@geno))){
        valid <- FALSE
        msg <- c(msg, "order of individuals differs in genotype and phenotype file!")
    }
    ## more snps in annotation, than genotyped
    if(length(unique(attr(object@geno,'anno')$snp))>ncol(object@geno)){
        valid <- FALSE
        msg <- c(msg, "annotation includes more SNPs than genotyped!")
    }
    ## SNPs in annotation file, that are not in genotype file (too big?)
    if(!all(unique(attr(dat,'anno')$snp) %in% colnames(object@geno))){
        valid <- FALSE
        msg <- c(msg, "there are SNPs in the annotation file that have no genotypes!")
    }
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

