#################################################
#
# GWASdata object functions
#
#################################################

#' @include pathway.r
NULL

#################################################

#' S4 class for an object representing a Genome-wide Assocaition Study.
#'
#' @slot geno An object of any type, including genotype information. The format
#' needs to be one line per individual and on colum per SNP in minor-allele 
#' coding (0,1,2). Other values between 0 and 2, as from impute dosages, are 
#' allowed. Missing values must be imputed prior to creation of a \code{GWASdata} object. 
#' @slot anno A \code{data.frame} mapping SNPs to genes and genes to
#' pathways. Needs to include the columns 'pathway' (pathway ID, e.g. hsa
#' number from KEGG database), 'gene' (gene name (hgnc_symbol)), 'chr'
#' (chromosome), 'snp' (rsnumber) and 'position' (base pair position of SNP).
#' @slot pheno A \code{data.frame} specifying individual IDs, phenotypes and
#' covariates to be included in the regression model e.g. ID, pheno, sex,
#' pack.years. Note: IDs have to be in the first column!
#' @slot desc A \code{character} giving the GWAS description, e.g. name of study.  
#' @examples
#' # create gwas data object
#' data(pheno)
#' data(geno)
#' data(anno)
#' gwas <- new('GWASdata', pheno=pheno, geno=geno, anno=anno, desc="some study") 
#' @author Juliane Manitz, Stefanie Friedrichs
#' @export GWASdata
#' @import methods
#' @rdname GWASdata-class
GWASdata <- setClass('GWASdata',
                     representation(geno ="ANY", anno="data.frame",
                                    pheno='data.frame', desc='character'))
    ## validy checks
    setValidity('GWASdata', function(object){
    msg   <- NULL
    valid <- TRUE
    ## check genotype has colnames (=SNP names) and rownames (individuals)
    if(any(is.null(colnames(object@geno)), is.null(rownames(object@geno)))){
        valid <- FALSE
        msg   <- c(msg, "genotypes need specified column- and/or rownames!")
    }
    ## check if genotypes include missings
     if(sum(is.na(object@geno))>0){
         valid <- FALSE
         msg   <- c(msg, "genotypes include missing values that need to be imputed!")
     }
    ## check phenotypes
    if(length(object@pheno)>0){
      ## phenotypes for more individuals than have genotypes
      if(length(unique(object@pheno[,1])) > length(row.names(object@geno))){
        valid <- FALSE
        msg <- c(msg,
                 "phenotypes exist for more individuals than have genotypes!")
      }
      ## check order of individuals in genotypes and phenotypes 
      if(!all(as.character(rownames(object@pheno)) %in% row.names(object@geno))){
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
        msg   <- c(msg, paste("column names in annotation are incorrect! They need to be:",anno_names))
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
#' \code{'GWASdata'} is a GWASdata object constructor.
#' @param geno An object of any type, including the genotype information.  
#' @param anno A \code{data.frame} containing the annotation file for the 
#' \code{GWASdata} object.
#' @param pheno A \code{data.frame} specifying individual IDs, phenotypes and
#' covariates to be included in the regression model.
#' @param desc A \code{character} giving the GWAS description, e.g. name of study. 
#' @param ... Further arguments can be added to the function.
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
#' read genotype data from file to one of several available objects, which 
#' can be passed to a GWASdata object \code{\link{GWASdata}}.
#'
#' @param file.path \code{character} giving the path to the data file to be read
#' @param save.path \code{character} containing the path for the backingfile
#' @param sep \code{character}. A field delimeter. See 
#' \code{\link[bigmemory]{read.big.matrix}} for details. 
#' @param header \code{logical}. Does the data set contain column names?
#' @param use.fread \code{logical}. Should the dataset be read using the function 
#' \code{\link[data.table]{fread}} \code{fread} from package \pkg{data.table}?
#' @param use.big \code{logical}. Should the dataset be read using the function 
#' \code{\link[bigmemory]{read.big.matrix}} from package \pkg{bigmemory}?
#' @param row.names \code{logical}. Does the dataset include rownames?
#' @param ... further arguments to be passed to \code{read_geno}.
#' @details If the data set contains rownames specified, set option \code{has.row.names = TRUE}.
#'
#' @importFrom bigmemory as.big.matrix read.big.matrix
#' @importFrom data.table fread
#' 
## @rdname read_geno
## @name read_geno
#' @aliases read_geno character
#' @examples 
#' \dontrun{
#' path <- system.file("extdata", "geno.txt", package = "kangar00")
#' geno <- read_geno(path, save.path = getwd(), sep = " ", use.fread = FALSE, row.names = FALSE)
#' }
#' @import tools
#' @export
setMethod("read_geno",
          signature="character",
                     #save.path="character", sep="character",
                     #header="logical", use.fread="logical"),
          definition = function(file.path, save.path = NULL, sep = " ",
                                header = TRUE, use.fread = TRUE, use.big = FALSE,
                                row.names = FALSE, ...) {

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
            file.name <- utils::tail(unlist(strsplit(file.path, "/")), 1)
            file.name <- unlist(strsplit(file.name, "[.]"))
            file.name <- paste(file.name[-length(file.name)], collapse=".")
            if(is.null(save.path)){
              save.path <- normalizePath(dirname(file.path))
              save.file <- paste(file.name, ".bin", sep="")
              warning(paste("Default save.path", save.path, "was created!", sep=" "))
            }else{
              save.file <- paste0(save.path, "/", file.name, ".bin", sep="")
            }
            
            ## Step 1.4: Check read options
            if(use.fread == TRUE & use.big == TRUE){
              stop(paste("You can not use the functions fread from package data.table and read.big.matrix from package bigmemory at the same time!"))
            }
            
            # Step 1.5: Check rownames and fread
            if(use.fread == TRUE & row.names == TRUE){
              stop(paste("The function fread from package data.table can not handle row.names! Please try to use the function read.big.matrix from package bigmemory or the read.table function!"))
            }

            # Step 2: Read in file according to format
            cat("Loading data. This might take a while depending on the size of the
                data set.")

            if (fileFormat == "gz"){
              stop(paste(fileFormat, "currently not supported", sep=" "))
              
              # Step 2.1
              cat("Reading in huge beagle files may fail due to memory limits.
                  If this is the case convert your beagle file in a .txt-file and try again. \n")
              if(use.fread){
                cat("Loading data via fread function from package data.table. If problems occur, your file will automatrically be read using the function read.table instead. \n")
                tryCatch({
                  gwasGeno <- data.table::fread(file.path, header = TRUE)
                }, warning = function(wr){
                  cat("fread caused a warning. Please check if your file has been read in correctly! \n")
                  print(wr)
                }, error = function(er){
                  cat("fread has stopped due to an error. \n")
                  print(er)
                  cat("Try to read your file using the read.table function. Attention: This function my be very slow on large files! \n")
                  tryCatch({
                    gwasGeno <- utils::read.table(file.path, header = TRUE)
                  }, warning = function(w){
                    cat("read.table caused a warning, Please check if your file has been read in correctly! \n")
                    print(w)
                  }, error = function(e){
                    cat("Also the read.table function has stopped due to an error. Try to convert your file in a .txt-file and try again. \n")
                    print(e)
                  }
                    )
                }
                  )
                gwasGeno <- bigmemory::as.big.matrix(gwasGeno)
              } else{
                gwasGeno <- utils::read.table(file.path, header=TRUE)
                gwasGeno <- bigmemory::as.big.matrix(gwasGeno)
              }
            } else if (fileFormat == "mldose"){
              stop(paste(fileFormat, "currently not supported!", sep=" "))
              cat("Reading in huge MACH files may fail due to memory limits.
                  If this is the case convert your MACH file in a .txt-file and try again.")
              if(use.fread){
                cat("Loading data using function fread. If this leads to problems set use.fread = FALSE")
                gwasGeno <- data.table::fread(file.path, header = TRUE)
                gwasGeno <- bigmemory::as.big.matrix(gwasGeno)

              } else{
                gwasGeno <- utils::read.table(file.path, header=TRUE)
                gwasGeno <- bigmemory::as.big.matrix(gwasGeno)
              }
            } else if (fileFormat == "impute2"){
              stop(paste(fileFormat, "currently not supported!", sep=" "))
              cat("Reading in huge IMPUTE2 files may fail due to memory limits.
                  If this is the case convert your IMPUTE2 file in a .txt-file and try again.")
              if(use.fread){
                cat("Loading data via fread. If this leads to problems set use.fread = FALSE")
                gwasGeno <- data.table::fread(file.path, header = TRUE)
                gwasGeno <- bigmemory::as.big.matrix(gwasGeno)

              } else{
                gwasGeno <- utils::read.table(file.path, header=TRUE)
                gwasGeno <- bigmemory::as.big.matrix(gwasGeno)
              }
            } else if (fileFormat == "txt"){
              if(use.big){
                options(bigmemory.allow.dimnames=TRUE)
                tryCatch({
                  print(save.file)
                  print(save.path)
                  if(row.names){
                    gwasGeno <- bigmemory::read.big.matrix(file.path, type='char',
                                                           backingfile = save.file,
                                                           backingpath = save.path,
                                                           descriptorfile = paste0(file.name, ".desc", sep=""),
                                                           sep = sep, header = header,
                                                           has.row.names = TRUE,...)
                  } else {
                    gwasGeno <- bigmemory::read.big.matrix(file.path, type='char',
                                                           backingfile = save.file,
                                                           backingpath = save.path,
                                                           descriptorfile = paste0(file.name, ".desc", sep=""),
                                                           sep = sep, header = header,
                                                           has.row.names = FALSE,...)
                  }
                  
                }, warning = function(w){
                  cat("read.big.matrix caused a warning! Please make sure your file has been read correctly!")
                  print(w)
                }, error = function(e){
                  cat("The function read.big.matrix has stopped due to an error!")
                  print(e)
                  cat("Try reading your file using the read.table function. Attention: This function is very slow for big files! \n")
                  tryCatch({
                    gwasGeno <- utils::read.table(file.path, header = TRUE)
                  }, warning = function(w){
                    cat("read.table caused a warning, Please make sure your file has been read in correctly! \n")
                    print(w)
                  }, error = function(e){
                    cat("Also the read.table function has stopped due to an error. \n")
                    print(e)
                  }
                  )
                }
                )
              } else if(use.fread){
                tryCatch({
                  gwasGeno <- data.table::fread(file.path, header = TRUE)
                }, warning = function(wr){
                  cat("The fread fucntion caused a warning. Please make sure your file has been read in correctly! \n")
                  print(wr)
                }, error = function(er){
                  cat("fread has stopped due to an error. \n")
                  print(er)
                  cat("Try reading in file with read.table. Attention: This function can be very slow for large files! \n")
                  tryCatch({
                    gwasGeno <- utils::read.table(file.path, header = TRUE)
                  }, warning = function(w){
                    cat("teh read.table function caused a warning, Please make sure your file has been read in correctly! \n")
                    print(w)
                  }, error = function(e){
                    cat("Also the read.table function has stopped due to an error. \n")
                    print(e)
                  }
                  )
              }
                )
              } else if(!use.fread & !use.big){
                cat("Try reading in file using the read.table function. Attention: This function may be very slow on big files! \n")
                gwasGeno <- utils::read.table(file.path, header = TRUE, sep = sep)
              }
            } else{
              stop("Unknown file format. Please only use
                   .txt-Files or output from MACH, Impute or Beagle!")
            }

            # Step 4: Check if output has right format
            # Step 4.1 Check if dataset contains data as expected
            # Step 4.1.1. Check if file has no row.names
            if(row.names){
              firstRow <- gwasGeno[ , 1]
              if(!is.character(stats::na.omit(firstRow)) || sum(firstRow > 2) == 0){
                warning("Your geno file doesn't seem to contain ID numbers. 
                Please make sure that the first row of your data contains ID 
                numbers according to your phenotype file! Otherwise use 
                row.names = FALSE!")
              }

              # Step 4.1.2 Check if the rest of the data is okay
              if(sum(stats::na.omit(gwasGeno[ , -1]) > 2) > 0){
                warning("Your genotypes contain values bigger than 2.")
              }
            } else{ # Step 4.2.1 Check if file has row.names
              possIDs <- rownames(gwasGeno)
              if(!is.character(stats::na.omit(possIDs)) || sum(possIDs > 2) == 0){
                warning("Your genotype file doesn't seem to contain ID numbers. 
                Please make sure that the first row of your data contains ID 
                numbers matching  your phenotype file!")
              }

              # Step 4.2.2 Check if rest of the data is okay
              if(sum(stats::na.omit(gwasGeno[,]) > 2) > 0){
                warning("Your genotype data seems to contain values bigger than 2.")
              }

            }

            # Step 5:
            return(gwasGeno)
})

#' \code{show} displays basic information on \code{\link{GWASdata}} object
#' @param object A \code{\link{GWASdata}} object.
#' @examples
#' # show and summary methods
#' gwas
## @author Juliane Manitz
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
              	print(utils::head(object@pheno))
              }
              invisible(NULL)
          })


## summary
setGeneric('summary', function(object, ...) standardGeneric('summary'))

#' \code{summary} summarizes the content of a \code{\link{GWASdata}} object 
#' and gives an overview about the information included in a 
#' \code{\link{GWASdata}} object. Summary statistics for phenotype and genotype 
#' data are calculated.
#'
#' @examples
## # summary method
## data(gwas) 
#' summary(gwas)
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

#' \code{GeneSNPsize} creates a \code{data.frame} of \code{\link{pathway}}
#' names with numbers of snps and genes in each \code{\link{pathway}}.
#'
#' @describeIn GWASdata creates a \code{data.frame} of \code{\link{pathway}} names with numbers 
#' of snps and genes in each pathway.
#' @export 
#' @aliases GeneSNPsize GWASdata
#' @examples
#' # SNPs and genes in pathway
## data(gwas) 
#' GeneSNPsize(gwas)
setMethod('GeneSNPsize', signature='GWASdata',
          definition <- function(object){
              nrsnps <- table(unique(object@anno[,c('pathway','gene','snp')])$pathway)
              nrgenes <- table(unique(object@anno[,c('pathway','gene')])$pathway)
              tab <- cbind(nrgenes,nrsnps)
	        colnames(tab) <- c('genes','SNPs')
          return(tab)
          })


#' An S4 class for an object assigning SNP positions to rs-numbers (for internal use) 
#'
#' @rdname snp_info-class
#' @slot info A \code{data.frame} including information on SNP positions
#' 
#' @author Stefanie Friedrichs
## @examples
## # compare with package data
## snp_info("rs10243170")
## data(rs10243170_info) 
#' @export snp_info
#' @import methods
#' @seealso \code{\link{pathway_info}}, \code{\link{get_anno}}
snp_info <- setClass('snp_info', slots=c(info='data.frame'))

setValidity('snp_info', function(object){  
	msg  <- NULL
	valid <- TRUE
 	if(!is.data.frame(object@info)){
	  valid=FALSE
	  msg <- c(msg, "the SNP_info object must include a data.frame")
	}
	if(!all.equal(colnames(object@info),c("chr","position","snp"))){
	  valid=FALSE
	  msg <- c(msg, "the included data.frame needs columns 'chr', 'position' 
    and 'snp'")
	}
})

setGeneric('snp_info', function(x, ...) standardGeneric('snp_info'))
#' This function gives for a \code{vector} of SNP identifiers the position of each SNP 
#' as extracted from the Ensemble database. The database is accessed via the
#' R-package \pkg{biomaRt}.
#'
#' @param x A \code{character} \code{vector} of SNP rsnumbers for which 
#' positions will be extracted.
#' @param ... further arguments can be added.
#' @return A \code{data.frame} including the SNP positions with columns
#' 'chromosome', 'position' and 'snp'. SNPs not found in the Ensemble database
#' will not be listed in the returned \code{snp_info} object, SNPs with multiple
#' positions  will appear several times.
#' @examples
## snp_info("rs234")
#' snp_info("rs10243170")
#'
## @author Stefanie Friedrichs
#' @import biomaRt
#' @rdname snp_info-class 
#' @export 
setMethod('snp_info', signature='character', 
          definition = function(x){        
  #set database; Homo sapiens Short Variation (SNPs and indels):
  # snp <- useMart("snp", dataset="hsapiens_snp") ##server unavailable!
  #listMarts(host = "jul2015.archive.ensembl.org")
  ensembl <- useMart(biomart="ENSEMBL_MART_SNP",dataset="hsapiens_snp", 
                     host = "jul2015.archive.ensembl.org")                         
  info <- getBM(attributes=c("chr_name","chrom_start","refsnp_id"),
                filters=c("snp_filter"),values=x, mart=ensembl)
  colnames(info) <- c("chr","position","snp")
  ret <- new('snp_info', info=info)
  return(ret)
})

#' \code{show} Shows basic information on \code{\link{snp_info}} object
#'
## @param object An \code{object} of class \code{\link{snp_info}}.
#' @return \code{show} Basic information on \code{\link{snp_info}} object.
## @examples
## # show
## data(rs10243170_info)
## rs10243170_info
## @author Stefanie Friedrichs
#' @export
#' @rdname snp_info-class
setMethod('show', signature='snp_info',
          definition = function(object){
            cat('An object of class ', class(object), '\n', sep='')
            cat('Number of SNPs:', nrow(object@info),'\n')
            if(nrow(object@info)>6){ cat('First six rows: \n')
            print(object@info[1:6,])}else{print(object@info)}  
})

setGeneric('summary', function(object, ...) standardGeneric('summary'))
#' \code{summary} Summarizes information on \code{\link{snp_info}} object
#'
#' @param object An \code{object} of class \code{\link{snp_info}}.
#' @return \code{summary} Summarized information on \code{\link{snp_info}} object.
#' @examples
#' summary(rs10243170_info)
## @author Stefanie Friedrichs
#' @export
#' @rdname snp_info-class
setMethod('summary', signature='snp_info',
          definition = function(object){
            cat('An object of class ', class(object), '\n', sep='')
            cat('Number of SNPs:', nrow(object@info),'\n')
            if(nrow(object@info)>6){ cat('First six rows: \n')
            print(object@info[1:6,])}else{print(object@info)}
               
})

setGeneric('get_anno', function(object1, object2, ...) standardGeneric('get_anno'))
#' Annotates SNPs via genes to pathways
#'
#' A function to create the annotation for a \code{\link{GWASdata}} object. 
#' It combines a \code{\link{snp_info}} and a \code{\link{pathway_info}}
#' object into an annotation \code{data.frame} used for \code{\link{pathway}} 
#' analysis on GWAS. SNPs are assigned to pathways via gene membership.
#'
#' @param object1 A \code{\link{snp_info}} object with SNP information as returned
#' by the \code{\link{snp_info}} function. The included \code{data frame} contains
#' the columns 'chr', 'position' and 'snp'.
#' @param object2 A \code{\link{pathway_info}} object with information on genes
#' contained in pathways. It is created by the \code{\link{pathway_info}} 
#' function and contains a \code{data frame} with columns 
#' 'pathway', 'gene_start', 'gene_end', 'chr', 'gene'.
#' @param ... further argdata(hsa04020)
#' @return A \code{data.frame} mapping SNPs to genes and genes to
#' pathways. It includes the columns 'pathway', 'gene', 'chr', 'snp' and 
#' 'position'.
#' 
#' @examples
#' data(hsa04022_info)  # pathway_info('hsa04020')
#' data(rs10243170_info)# snp_info("rs10243170")
#' get_anno(rs10243170_info, hsa04022_info)
#'
#' @author Stefanie Friedrichs, Saskia Freytag, Ngoc-Thuy Ha
#' @seealso \code{\link{snp_info}}, \code{\link{pathway_info}}
## @rdname snp_info-methods
#' @aliases get_anno
#' @import sqldf
#' @export
setMethod('get_anno', signature=c('snp_info','pathway_info'),
          definition <- function(object1, object2, ...) {
          
  if (!inherits(object1, "snp_info"))
      stop("First argument needs to be a snp_info object!")
  if (!inherits(object2, "pathway_info"))  
      stop("Second argument has to be a pathway_info object!")
  if (!(all.equal(colnames(object1@info),c("chr", "position", "snp"))))
      stop("The data frame included in the snp_info object has incorrect column names! 
            Column names need to be: 'chr', 'positon' and 'snp'!")
  if (!(all.equal(colnames(object2@info),c("pathway","gene_start",
      "gene_end","chr","gene"))))
      stop("The data frame included in the pathway_info object has incorrect column names! Column names need to be: 'pathway', 'gene_start', 'gene_end', 'chr' and 'gene'!")  
                              
  object1@info$chr <- as.factor(object1@info$chr)
  object2@info$chr <- as.factor(object2@info$chr)  
  snp_info     <- split(object1@info, object1@info$chr)
  pathway_info <- split(object2@info, object2@info$chr) 
  
  list_out <- list()
  for(i in 1:length(names(pathway_info))){
    x <- pathway_info[[i]] 
    y <- try(snp_info[[which(names(snp_info)==names(pathway_info)[i])]], silent=T)
    if(is.null(dim(y))) next
    list_out[[i]] <- sqldf("select x.pathway, x.gene,
                                   y.chr ,y.snp, y.position
                                   from x x, y y
                                   where x.gene_start<=y.position AND
                                   x.gene_end>=y.position")
    remove(x,y)
  }      
  anno <- do.call("rbind", lapply(list_out, data.frame, stringsAsFactors=FALSE))
  return(anno)
})
