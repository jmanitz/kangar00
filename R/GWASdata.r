#################################################
#
# GWASdata object functions
#
#################################################

#' An S4 class defining an object to represent a Genome-wide Assocaition Study.
#'
#' @slot geno An \code{big.matrix} object including genotype information.
#' @slot anno A \code{data.frame} mapping SNPs to genes and genes to
#' pathways. Needs to include the columns 'pathway' (pathway ID, e.g. hsa
#' number from KEGG database), 'gene' (gene name (hgnc_symbol)), 'chr'
#' (chromosome), 'snp' (rsnumber) and 'position' (base pair position of SNP).
#' @slot pheno A \code{data.frame} specifying individual IDs, phenotypes and
#' covariates to be included in the regression model e.g. ID, pheno, sex,
#' pack.years. Note: IDs have to be in the first column!
#' @slot desc A \code{character} giving the GWAS description, e.g. name of study
#' @examples
#' data(pheno)
#' data(geno)
#' #gwas <- new('GWASdata', pheno=pheno, geno=geno, desc="some study") ### ERROR
#' @author Juliane Manitz, Stefanie Friedrichs
#' @exportClass GWASdata
#' @export GWASdata
#' @import methods
GWASdata <- setClass('GWASdata', slots=c(geno="big.matrix", anno = "data.frame",
	                                 pheno='data.frame', desc='character'))
    ## validy checks
    setValidity('GWASdata', function(object){
    msg   <- NULL
    valid <- TRUE
    ## check genotype has colnames (=SNP names) and rownames (individuals)
    if(any(is.null(colnames(object@geno)), is.null(rownames(object@geno)))){
        valid <- FALSE
        msg   <- c(msg, "object geno needs specified col- and/or rownames")
    }
    ## check whether GWASdata@geno has missings
    # <FIXME> Currently not working
    # if(sum(is.na(object@geno))>0){
    #     valid <- FALSE
    #     msg   <- c(msg, "genotypes include missing values and need to be imputed!")
    # }
    ## check phenotypes
    if(length(object@pheno)>0){
      ## phenotypes for more individuals than have genotypes
      if(length(unique(object@pheno[,1])) > length(row.names(object@geno))){
        valid <- FALSE
        msg <- c(msg,
                 "phenotypes exist for more individuals than have genotypes!")
      }
      ## check order of individuals in genotypes and phenotypes
      if(!all.equal(as.character(object@pheno[,1]),row.names(object@geno))){
        valid <- FALSE
        msg <- c(msg,
                 "order of individuals differs in genotype and phenotype file!")
      }
    }else{
      message("Note that phenotypes are not specified.")
    }
    ## check annotation file: variable names
    anno_names <- c('pathway','gene','chr','snp','position')
    if(!all(anno_names %in% colnames(object@anno))){
        valid <- FALSE
        msg   <- c(msg, paste("anno variable names do not match:",anno_names))
    }
    ## more snps in annotation, than genotyped: <FIXME> Why do we needthis? <FIXME>
    if(length(unique(object@anno$snp)) > ncol(object@geno)){
        valid <- FALSE
        msg <- c(msg, "annotation includes more SNPs than genotyped!")
    }
    ## SNPs in annotation file, that are not in genotype file (too big?)
    if(!all(unique(object@anno$snp) %in% colnames(object@geno))){
        valid <- FALSE
        msg <- c(msg, "there are SNPs in the annotation file that have no genotypes!")
    }
    if(valid) TRUE else msg
})

#' \code{'GWASdata'} is a GWASdata object constructor
#'
#' @param geno an \code{ff} data frame including genotype information.
#' @param anno a \code{data.frame} mapping SNPs to genes and genes to
#' pathways. Needs to include the columns 'pathway' (pathway ID, e.g. hsa
#' number from KEGG database), 'gene' (gene name (hgnc_symbol)), 'chr'
#' (chromosome), 'snp' (rsnumber) and 'position' (base pair position of SNP).
#' @param pheno A \code{data.frame} specifying individual IDs, phenotypes and
#' covariates to be included in the regression model e.g. ID, pheno, sex,
#' pack.years. Note: IDs have to be in the first column. Default is \code{NULL}.
#' @param desc A \code{character} giving the GWAS description, e.g. name of study. Defaul is ''.
#'
#' #@author Juliane Manitz
#' @export
#' @rdname GWASdata
#' @aliases show,GWASdata,ANY-method
setMethod('GWASdata',
       definition = function(geno, anno, pheno = NULL, desc = ''){
       ## create GWASdata object
       new('GWASdata', geno = geno, anno = anno, pheno = as.data.frame(pheno),
           desc = desc)
})
# GWASdata object constructor
setGeneric('GWASdata')#, function(object, ...) standardGeneric('GWASdata'))

# read genotype data from file
setGeneric('read_geno', function(x, ...) standardGeneric('read_geno'))
#' read genotype data from file to bigmemory object, which can be passed to a  GWASdata object
#'
#' @param file.path character, which contains the path to the data file to be read
#' @param save.path character, which contains the path for the backingfile
#' @param sep character. A field delimeter. See \code{bigmemory::read.big.matrix} for
#'  details.
#' @param header logical. Does the data set contain column names?
#' @param ... further arguments to be passed to \code{bigmemory::read.big.matrix}.
#' @details If the data set contains rownames specified, set option \code{has.row.names = TRUE}.
#'
#' @rdname read_geno
#' @seealso \code{\link{GWASdata-class}}
#' @import bigmemory
#' @export
setMethod('read_geno', signature='character',
       definition = function(x, save.path = NULL, sep = " ",
                             header = TRUE, ...) {
       cat("Loading data. This might take a while, depending on the size of the data set.\n")
            
       ## backing file path
       if(is.null(save.path))
          save.path <- paste0(x, '.bin')

       ## read data frame
       options(bigmemory.allow.dimnames=TRUE)

## <FIXME>:
## we should explain better what is expected (i.e., what data files can be
## imported) and show all possible options (type, sep, header, ... to the user)
## and not only implicitely via ... .
           df <- read.big.matrix(x, type='char',
                                 backingfile = save.path,
                                 descriptorfile = paste0(save.path, '.desc'),
                                 sep = sep, header = header, ...)
           rownames(df) <- paste0("ind", 1:nrow(df))
           return(df)
})

#' \code{show} Shows basic information on \code{GWASdata} object
#'
#' @param object A \code{GWASdata} object.
#' @return \code{show} This function shows the phenotype information and the object description included in a \code{\link{GWASdata}} object.
#'
#' @examples
#' #data(gwas) <FIXME> define example with new structure
#' #show(gwas)
#'
#' #@author Juliane Manitz
#' @export
#' @rdname GWASdata-class
#' @aliases show,GWASdata,ANY-method
setMethod('show', signature='GWASdata',
          definition = function(object){
              cat('An object of class ', class(object),
                  ' from ',object@desc,'\n\n',sep='')

              if(length(object@pheno) == 0){
		cat('No phenotypes specified.\n')
              }else{ # summary of phenotype and covariate data
                cat('Phenotypes for ', dim(object@pheno)[1],
                ' individuals: \n',sep='')
              	print(head(object@pheno))
              }
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
#' #data(gwas)   <FIXME> define example with new structure
#' #summary(gwas)
#'
#' #@author Juliane Manitz
#' @export
#' @rdname GWASdata-class
#' @aliases summary,GWASdata,ANY-method
setMethod('summary', signature='GWASdata',
          definition = function(object){
              cat('An object of class ', class(object), ' from ',object@desc,'\n\n',sep='')
              if(length(object@pheno) == 0){
		cat('No phenotypes specified.\n')
              }else{ # summary of phenotype and covariate data
	        cat('Phenotypes for ', dim(object@pheno)[1],
                    ' individuals: \n\n',sep='')
                print(summary(object@pheno))
              }
              cat('\nGenotypes: \n\n')
	      cat('Total number of genes and SNPs per pathway:\n')
              print(GeneSNPsize(object)) # overview of pathway - gene - snp mapping
# 	      cat('\n')
#	      print(object@geno)
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
#' #data(gwas) <FIXME> define example with new structure
#' #GeneSNPsize(gwas)
#'
#' #@author Juliane Manitz
#' @export
#' @rdname GWASdata-class
#' @aliases GeneSNPsize,GWASdata,ANY-method
setMethod('GeneSNPsize', signature='GWASdata',
          definition <- function(object){
              nrsnps <- table(unique(object@anno[,c('pathway','gene','snp')])$pathway)
              nrgenes <- table(unique(object@anno[,c('pathway','gene')])$pathway)
              tab <- cbind(nrgenes,nrsnps)
	        colnames(tab) <- c('genes','SNPs')
          return(tab)
          })


# <FIXME> create object snp_info and corresponding constructor <FIXME>
#setGeneric('snp_info', function(x, ...) standardGeneric('snp_info'))
#' Get SNP positions
#'
#' This function gives for a vector of SNPs the position of each SNP as
#' extracted from the Ensemble database. The database is accessed via the
#' R-package \code{biomaRt}.
#'
#' @param x a vector of SNP rsnumbers for which positions will be extracted. <FIXME> what type of vector? specify signature
#' @param ... further arguments can be added.
#' @return a \code{data.frame} including the SNP positions with columns
#' chromosome, position and rsnumber. SNPs not found in the Ensemble database
#' not be listed in the returned \code{data.frame}, SNPs with multiple positions
#' will appear several times.
#' @examples
#' # snp_info("rs234") <FIXME>
#'
#' @author Stefanie Friedrichs
#' @import biomaRt
setMethod('snp_info', signature='character',
          definition <- function(x, ...) {
  #set database; Homo sapiens Short Variation (SNPs and indels):
  snp <- useMart("snp", dataset="hsapiens_snp")
  snp_info <- getBM(attributes=c("chr_name","chrom_start","refsnp_id"),
              filters=c("snp_filter"),values=x, mart=snp)
  colnames(snp_info) <- c("chr","position","rsnumber")
  return(snp_info)
}

setGeneric('get_anno', function(obecjt, ...) standardGeneric('get_anno'))
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
#' #### missing example #### <FIXME>
#'
#' @author Stefanie Friedrichs, Saskia Freytag, Ngoc-Thuy Ha
#' @import sqldf
#setMethod('get_anno', signature='snp_info',
#          definition <- function(object, ...) {
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


