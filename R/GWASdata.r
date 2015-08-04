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

#' An S4 class defining an object to represent a Genome-wide Assocaition Study.
#'
#' @slot pheno A \code{data.frame} specifying individual IDs, phenotypes and
#' covariates to be included in the regression model e.g. ID, pheno, sex,
#' pack.years. Note: IDs have to be in the first column!
#' @slot geno An \code{ffdf} data frame including genotype information.
#' Has an attribute anno which is a \code{data.frame} mapping SNPs to genes and genes to
#' pathways. Needs to include the columns 'pathway' (pathway ID, e.g. hsa
#' number from KEGG database), 'gene' (gene name (hgnc_symbol)), 'chr'
#' (chromosome), 'snp' (rsnumber) and 'position' (base pair position of SNP).
#' @slot desc A \code{character} giving the GWAS description, e.g. name of study
#' @examples
#' data(pheno)
#' data(geno)
#' #gwas <- new('GWASdata', pheno=pheno, geno=geno, desc="some study") ### ERROR
#' @author Juliane Manitz, Stefanie Friedrichs
#' @exportClass GWASdata
#' @export GWASdata
#' @import methods
GWASdata <- setClass('GWASdata', slots=c(pheno='data.frame', geno='ffdf', desc='character'))
    ## validy checks
    setValidity('GWASdata', function(object){
    msg   <- NULL
    valid <- TRUE
    ## check annotation file
    if(is.null(attr(object@geno,"anno"))){
        valid <- FALSE
        msg   <- c(msg, "ffdf object geno needs an additional attribute:
                 data.frame 'anno'")
    }
    if(!is.data.frame(attr(object@geno,"anno"))){
  	    valid <- FALSE
        msg   <- c(msg, "geno attribute anno needs to be a data frame")
    }
    ## check anno data frame variable names
    anno_names <- c('pathway','gene','chr','snp','position')
    if(!all(anno_names %in% colnames(attr(object@geno,'anno')))){
        valid <- FALSE
        msg   <- c(msg, paste("anno variable names do not match:",anno_names))
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


#' \code{show} Shows basic information on \code{GWASdata} object
#'
#' @param object A \code{GWASdata} object.
#' @return \code{show} This function shows the phenotype information and the object description included in a \code{\link{GWASdata}} object.
#'
#' @examples
#' data(gwas)
#' show(gwas)
#'
#' #@author Juliane Manitz
#' @export
#' @rdname GWASdata-class
#' @aliases show,GWASdata,ANY-method
setMethod('show', signature='GWASdata',
          definition = function(object){
              cat('An object of class ', class(object), ' from ',object@desc,'\n\n',sep='')
              cat('Phenotypes for ', dim(object@pheno)[1], ' individuals: \n',sep='')
              print(head(object@pheno))
              invisible(NULL)
          })


## summary
setGeneric('summary', function(object, ...) standardGeneric('summary'))

#' \code{summary} Summarizes the content of a \code{GWASdata} object
#'
#' @param object A \code{GWASdata} object.
#' @return \code{summary} This function gives an overview about the information included in a \code{\link{GWASdata}} object. Summary statistics for phenotype and genotype data are calculated.
#'
#' @examples
#' # data(gwas)   #### ERROR
#' # summary(gwas) ### ERROR
#'
#' #@author Juliane Manitz
#' @export
#' @rdname GWASdata-class
#' @aliases summary,GWASdata,ANY-method
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

## GeneSNPsize
setGeneric('GeneSNPsize', function(object, ...) standardGeneric('GeneSNPsize'))

#' \code{GeneSNPsize} Counts the number of Genes and SPNs in each pathway
#'
#' @param object A \code{GWASdata} object.
#' @return \code{GeneSNPsize} Creates a list of pathway names with numbers of snps and genes in the pathway.
#'
#' @examples
#' #data(gwas)### ERROR
#' #GeneSNPsize(gwas) ### ERROR
#'
#' #@author Juliane Manitz
#' @export
#' @rdname GWASdata-class
#' @aliases GeneSNPsize,GWASdata,ANY-method
setMethod('GeneSNPsize', signature='GWASdata',
          definition <- function(object){
              anno <- attr(object@geno, 'anno')
              nrsnps <- table(unique(anno[,c('pathway','gene','snp')])$pathway)
              nrgenes <- table(unique(anno[,c('pathway','gene')])$pathway)
              tab <- cbind(nrgenes,nrsnps)
	        colnames(tab) <- c('genes','SNPs')
          return(tab)
          })


#' Get SNP positions
#'
#' This function gives for a vector of SNPs the position of each SNP as
#' extracted from the Ensemble database. The database is accessed via the
#' R-package \code{biomaRt}.
#'
#' @param snps A vector of SNP rsnumbers for which positions will be extracted.
#' @param ... further arguments can be added.
#' @return A \code{data.frame} including the SNP positions with columns
#' chromosome, position and rsnumber. SNPs not found in the Ensemble database
#' not be listed in the returned \code{data.frame}, SNPs with multiple positions
#' will appear several times.
#' @examples
#' # snp_info("rs234") ### ERROR
#'
#' @author Stefanie Friedrichs
#' @import biomaRt
snp_info <- function(snps, ...) {
  #set database; Homo sapiens Short Variation (SNPs and indels):
  snp <- useMart("snp", dataset="hsapiens_snp")
  snp_info <- getBM(attributes=c("chr_name","chrom_start","refsnp_id"),
              filters=c("snp_filter"),values=snps, mart=snp)
  colnames(snp_info) <- c("chr","position","rsnumber")
  return(snp_info)
}

#' Create annotation for GWASdata object
#'
#' This function gives for a vector of SNPs the position of each SNP as
#' extracted from the Ensemble database. The database is accessed via the
#' R-package \code{biomaRt}.
#'
#' @param snp_info A \code{data frame} with SNP information as returned by
#' \code{\link{snp_info}}. The \code{data frame} has to contain columns "chr",
#' "position" and "rsnumber".
#' @param pathway_info A \code{data frame} with information on the genes forming
#' the pathway. Output from \code{\link{pathway_info}}.
#' @param ... further arguments can be added.
#' @return A \code{data.frame} mapping SNPs to genes and genes to pathways. Includes
#' the columns "pathway", "gene", "chr", "snp" and "position".
#' @examples
#' #### missing example ####
#'
#' @author Stefanie Friedrichs, Saskia Freytag, Ngoc-Thuy Ha
#' @import sqldf
get_anno <- function(snp_info, pathway_info, ...){
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


