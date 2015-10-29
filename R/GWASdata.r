#################################################
#
# GWASdata object functions
#
#################################################

#' An S4 class defining an object to represent a Genome-wide Assocaition Study.
#'
#' @slot geno an \code{big.matrix} object including genotype information.
#' @slot anno A \code{data.frame} mapping SNPs to genes and genes to
#' pathways. Needs to include the columns 'pathway' (pathway ID, e.g. hsa
#' number from KEGG database), 'gene' (gene name (hgnc_symbol)), 'chr'
#' (chromosome), 'snp' (rsnumber) and 'position' (base pair position of SNP).
#' @slot pheno A \code{data.frame} specifying individual IDs, phenotypes and
#' covariates to be included in the regression model e.g. ID, pheno, sex,
#' pack.years. Note: IDs have to be in the first column!
#' @slot desc A \code{character} giving the GWAS description, e.g. name of study.  Default is ''.
#' @examples
#' data(pheno)
#' data(geno)
#' #gwas <- new('GWASdata', pheno=pheno, geno=geno, desc="some study") ### ERROR
#' @author Juliane Manitz, Stefanie Friedrichs
#' @export GWASdata
#' @import methods
#' @rdname GWASdata-class
GWASdata <- setClass('GWASdata', 
                     representation(geno ="big.matrix", anno="data.frame",  
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

# GWASdata object constructor
setGeneric('GWASdata', function(object, ...) standardGeneric('GWASdata'))
#' \code{'GWASdata'} is a GWASdata object constructor
#'
#' @export
#' @rdname GWASdata-class
setMethod('GWASdata',
       definition = function(geno, anno, pheno = NULL, desc = ''){
       ## create GWASdata object
       new('GWASdata', geno = geno, anno = anno, pheno = as.data.frame(pheno),
           desc = desc)
})

# read genotype data from file
setGeneric('read_geno', function(file.path, ...) standardGeneric('read_geno'))
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
#' @import bigmemory, tools
#' @export
setMethod("read_geno", 
          signature="character", 
                     #save.path="character", sep="character",
                     #header="logical", use.fread="logical"),
          definition = function(file.path, save.path = NULL, sep = " ",
                                header = TRUE, use.fread = TRUE,
                                row.names = TRUE, ...) {
            
            # Step 1: Check if file meets requirements
            # Step 1.1: Check if file exists
            if(!file.exists(file.path)){
              stop(paste("File", file.path, "does not exist!", sep=" "))
            }
            
            # Step 1.2: Check for right input file format
            fileFormat <- unlist(strsplit(file.path, "[.]"))
            fileFormat <- fileFormat[length(fileFormat)]
            
            acceptedFormats <- c("txt", "impute2", "gz", "mldose")
            if(!is.element(fileFormat, acceptedFormats)){
              stop(paste(fileFormat, "as file format is not accepted!", sep=" "))
            }
            
            ## Step 1.3: Check for backing file path otherwise use default
            file.name <- tail(unlist(strsplit(file.path, "/")), 1)
            file.name <- unlist(strsplit(file.name, "[.]"))
            file.name <- paste(file.name[-length(file.name)], collapse=".")
            if(is.null(save.path)){
              save.path <- normalizePath(dirname(file.path))
              save.file <- paste(file.name, ".bin", sep="")
              warning(paste("Default save.path", save.path, "was created!", sep=" "))
            }else{
              save.file <- paste0(save.path, "/", file.name, ".bin", sep="")
            }
            
            # Step 2: Read in file according to format
            cat("Loading data. This might take a while depending on the size of the 
                data set.")
            
            if (fileFormat == "gz"){
              # Step 2.1 
              car("Reading in huge beagle files may fail due to memory limits. 
                  If this is the case convert your beagle file in a .txt-file and try again. \n")
              if(use.fread){
                car("Loading data via fread. If this leads to problems set the function will try
                    automatically to load file via read.table. \n")
                tryCatch({
                  gwasGeno <- fread(file.path, header = TRUE)
                }, warning = function(wr){
                  cat("fread caused a warning. Please make sure your file has been read in correctly! \n")
                  print(wr)
                }, error = function(er){
                  cat("fread has stopped due to an error. \n")
                  print(er)
                  cat("Try reading in file with read.table. Attention: This function is very slow! \n")
                  tryCatch({
                    gwasGeno <- read.table(file.path, header = TRUE)
                  }, warning = function(w){
                    cat("read.table caused a warning, Please make sure your file has been read in correctly! \n")
                    print(w)
                  }, error = function(e){
                    cat("Also read.table has stopped due to an error. Try to convert your file in a .txt-file
                        and try again. \n")
                    print(e)
                  }
                    )
                }
                  )
                gwasGeno <- as.big.matrix(gwasGeno)
              } else{
                gwasGeno <- read.table(file.path, header=TRUE)
                gwasGeno <- as.big.matrix(gwasGeno)
              }
            } else if (fileFormat == "mldose"){
              cat("Reading in huge MACH files may fail due to memory limits. 
                  If this is the case convert your MACH file in a .txt-file and try again.")
              if(use.fread){
                cat("Loading data via fread. If this leads to problems set use.fread = FALSE")
                gwasGeno <- fread(file.path, header = TRUE)
                gwasGeno <- as.big.matrix(gwasGeno)
                
              } else{
                gwasGeno <- read.table(file.path, header=TRUE)
                gwasGeno <- as.big.matrix(gwasGeno)
              }
            } else if (fileFormat == "impute2"){
              cat("Reading in huge IMPUTE2 files may fail due to memory limits. 
                  If this is the case convert your IMPUTE2 file in a .txt-file and try again.")
              if(use.fread){
                cat("Loading data via fread. If this leads to problems set use.fread = FALSE")
                gwasGeno <- fread(file.path, header = TRUE)
                gwasGeno <- as.big.matrix(gwasGeno)
                
              } else{
                gwasGeno <- read.table(file.path, header=TRUE)
                gwasGeno <- as.big.matrix(gwasGeno)
              }
            } else if (fileFormat == "txt"){
              options(bigmemory.allow.dimnames=TRUE)
              tryCatch({
                print(save.file)
                print(save.path)
                if(row.names){
                  gwasGeno <- read.big.matrix(file.path, type='char',
                                              backingfile = save.file,
                                              backingpath = save.path,
                                              descriptorfile = paste0(file.name, ".desc", sep=""),
                                              sep = sep, header = header,
                                              has.row.names = TRUE,...)
                } else {
                  gwasGeno <- read.big.matrix(file.path, type='char',
                                              backingfile = save.file,
                                              backingpath = save.path,
                                              descriptorfile = paste0(file.name, ".desc", sep=""),
                                              sep = sep, header = header,
                                              has.row.names = FALSE,...)
                }
             
              }, warning = function(w){
                cat("read.big.matrix caused a warning!")
                print(w)
              }, error = function(e){
                cat("read.big.matrix has stopped due to an error!")
                print(e)
              }
                )
              
            } else{
              stop("Unknown file format. Please only use 
                   .txt-Files or output from MACH, Impute or Beagle!")
            }
            
            # Step 4: Check if output has right format
            # Step 4.1 Check if dataset contains data as expected
            # Step 4.1.1. Check if file has no row.names
            if(!row.names){
              firstRow <- gwasGeno[ , 1]
              if(!is.character(na.omit(firstRow)) || sum(firstRow > 2) == 0){
                warning("Your geno file doesn't seem to contain ID numbers. Please make sure that
                   the first row of your data contains ID numbers according to your phenotype file!")
              }
              
              # Step 4.1.2 Check if the rest of the data is okay
              if(sum(na.omit(gwasGeno[ , -1]) > 2) > 0){
                warning("Your geno data seems to contain values bigger than 2.")
              }
            } else{ # Step 4.2.1 Check if file has row.names
              possIDs <- rownames(gwasGeno)
              if(!is.character(na.omit(possIDs)) || sum(possIDs > 2) == 0){
                warning("Your geno file doesn't seem to contain ID numbers. Please make sure that
                   the first row of your data contains ID numbers according to your phenotype file!")
              }
              
              # Step 4.2.2 Check if rest of the data is okay
              if(sum(na.omit(as.matrix(gwasGeno)) > 2) > 0){
                warning("Your geno data seems to contain values bigger than 2.")
              }
              
            }
            
            # Step 5: Change geno object to big.matrix.object
            return(gwasGeno)
})

#' \code{show} displays basic information on \code{GWASdata} object
#' @param object A \code{GWASdata} object.
#' #@author Juliane Manitz
#' @export
#' @rdname GWASdata-class
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

#' \code{summary} summarizes the content of a \code{GWASdata} object and gives an overview about the information included in a \code{\link{GWASdata}} object. Summary statistics for phenotype and genotype data are calculated.
#'
#' @examples
#' #data(gwas)   <FIXME> define example with new structure
#' #summary(gwas)
#'
#' @export
#' @rdname GWASdata-class
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

#' \code{GeneSNPsize} creates a \code{data.frame} of pathway names with numbers of snps and genes in each pathway.
#'
#' @examples
#' #data(gwas) <FIXME> define example with new structure
#' #GeneSNPsize(gwas)
#'
#' @export
#' @rdname GWASdata-class
setMethod('GeneSNPsize', signature='GWASdata',
          definition <- function(object){
              nrsnps <- table(unique(object@anno[,c('pathway','gene','snp')])$pathway)
              nrgenes <- table(unique(object@anno[,c('pathway','gene')])$pathway)
              tab <- cbind(nrgenes,nrsnps)
	        colnames(tab) <- c('genes','SNPs')
          return(tab)
          })


# <FIXME> create object snp_info and corresponding constructor <FIXME>
setGeneric('snp_info', function(x, ...) standardGeneric('snp_info'))
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
})

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

