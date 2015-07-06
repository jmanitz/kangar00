#################################################
#
# GWASdata object functions
#
#################################################

## make sure that class "ffdf" is formally defined
#' "ffdf" class for memory-efficient storage of genotype data frame on disk and fast access
#'
#' @name ffdf-class
#' @aliases ffdf
#' @family ffdf
#'
#' @import ff
#' @exportClass ffdf
#' @seealso ff::ffdf ff::read.table.ff
setOldClass('ffdf')

# object constructor
GWASdata <- setClass('GWASdata', slots=c(pheno='data.frame', geno='ffdf', desc='character'))

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
    if(length(unique(object@pheno[,1]))>length(row.names(object@geno))){
        valid <- FALSE
        msg <- c(msg, "phenotypes exist for more individuals than have genotypes!")
    }
    ## check order of individuals in genotypes and phenotypes
    if(!all.equal(as.character(object@pheno[,1]),row.names(object@geno))){
        valid <- FALSE
        msg <- c(msg, "order of individuals differs in genotype and phenotype file!")
    }
    ## more snps in annotation, than genotyped
    if(length(unique(attr(object@geno,'anno')$snp))>ncol(object@geno)){
        valid <- FALSE
        msg <- c(msg, "annotation includes more SNPs than genotyped!")
    }
    ## SNPs in annotation file, that are not in genotype file (too big?)
    if(!all(unique(attr(object@geno,'anno')$snp) %in% colnames(object@geno))){
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
              nrsnps <- table(unique(anno[,c('pathway','gene','snp')])$pathway)
              nrgenes <- table(unique(anno[,c('pathway','gene')])$pathway)
              tab <- cbind(nrgenes,nrsnps)
	      colnames(tab) <- c('genes','SNPs')
              return(tab)
          })

## Update/get snp position
snp_info <- function(snps, ...) {
## snp_info: vector of SNP rsnumbers for which positions will be extracted 
  #if (!inherits(snps, "vector"))
  #    stop("SNPs object is not a vector")  ?

  #set database; Homo sapiens Short Variation (SNPs and indels):
  snp <- useMart("snp", dataset="hsapiens_snp")
  #SNPs not found are not listed
  snp_info <- getBM(attributes=c("chr_name","chrom_start","refsnp_id"),
              filters=c("snp_filter"),values=snps, mart=snp)
  colnames(snp_info) <- c("chr","position","rsnumber")          
  return(snp_info) #data frame with chromosome, position and rsnumber
}

# create annotation for GWASdata object
get_anno <- function(snp_info, pathway_info, ...){
## snp_info : data frame with SNP information columns called "chr", 
## "position" and "rsnumber" are assumed
## pathway_info: dataframe, pathway-gene information file from getPathway
  if (!inherits(snp_info, "data.frame"))
      stop("SNP object is not a data frame")
  if (!inherits(pathway_info, "data.frame"))
      stop("pathway gene information object is not a data frame")
  if ( !(columns(snp_info)%in% c("chr", "position", "rsnumber")) )
      stop("SNP data frame needs columns for chromosome, positon and rsnumber")
  if ( !(columns(pathway_info)%in% c("pathway","gene_start",
                                    "gene_end","chr","gene")) )
      stop("pathway information data frame needs columns for pathway, 
            gene_start, gene_end, chr and gene")      
  pathwayanno$Chr <- as.factor(pathway_info$chr)
  snp_info$chr <- as.factor(snp_info$chr)
  pathway_info <- split(pathway_info, pathway_info$chr)
  snp_info <- split(snp_info, snp_info$chr)

  list_out <- list()
  for(i in 1:length(names(pathway_info))){
    x <- pathway_info[[i]]
    y <- try(snp_info[[which(names(snp_info)==names(pathway_info)[i])]], silent=T)
    if(is.null(dim(y))) next
    list_out[[i]] <- sqldf("select x.pathway, x.gene,
                                   y.chr ,y.rsnumber, y.position
                                   from x x, y y
                                   where x.gene_start>=y.position AND
                                   x.gene_end<=y.position")
    remove(x,y)
  }
  anno <- lapply(list_out,"names<-",
          value=c("pathway", "gene", "chr","snp", "position"))
  anno <- do.call("rbind", lapply(anno, data.frame, stringsAsFactors = FALSE))
  return(anno)
}


