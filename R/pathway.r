#################################################
#
# pathway object functions
#
#################################################

#' An S4 class to represent a gene-gene interaction network
#'
#' @rdname pathway-class
#'
#' @slot id A \code{character} repesenting the \code{\link{pathway}} id, 
#' e.g. hsa00100 as used in the KEGG database.
#' @slot adj A \code{matrix} respresenting the network adjacency matrix of dimension 
#' equaling the number of genes (1 interaction, 0 otherwise)
#' @slot sign A \code{numeric} \code{vector} indicating the interaction type for 
#' each link (1 activation, -1 inhibition) in the interaction network for the 
#' \code{\link{pathway}}.
#' 
#' @author Juliane Manitz, Stefanie Friedrichs, Patricia Burger
#'
#' @export
#' @import methods
pathway <- setClass('pathway',
                    slots=c(id='character', adj='matrix', sign='vector'))

# validy checks
setValidity('pathway', function(object){ 
	msg  <- NULL
	valid <- TRUE
	# adjacency includes only 1 and 0
	if(!all( object@adj %in% c(0,1))){
	  valid=FALSE
	  msg <- c(msg, "the adjacency matrix is not allowed to include other values 
           that zero and one")
	}
	# length(sign) = number of links
	if(length(object@sign)!=sum(object@adj!=0)){
	  valid <- FALSE
          msg <- c(msg, "the length of the sign vector does not correspond to 
                 the number of links")
        }
	# quadratic?
	if(nrow(object@adj)!=ncol(object@adj)){
	  valid <- FALSE
          msg <- c(msg, "adjacency matrix has to be quadratic")
        }
	# isSymmetric(adj)
	if(!isSymmetric(object@adj)){
	  valid <- FALSE
          msg <- c(msg, "adjacency matrix has to be symmetric")
        }
	if(valid) TRUE else msg
})

# pathway object constructor
setGeneric('pathway', function(object, ...) standardGeneric('pathway'))
#' \code{pathway} is the \code{\link{pathway}} object constructor.
#'
#' @param id A \code{character} repesenting the \code{\link{pathway}} id.  
#' @param adj A \code{matrix} respresenting the network adjacency matrix of dimension 
#' equaling the number of genes (1 interaction, 0 otherwise)
#' @param sign A \code{numeric} \code{vector} indicating the interaction type for  
#' each link (1 activation, -1 inhibition) in the interaction network for the \code{\link{pathway}}.
## @param ... Further arguments can be added to the function.
#' @export
#'
#' @examples
#' # pathway object constructor
#' pathway(id="hsa04022")
#' 
#' @rdname pathway-class
setMethod('pathway',
       definition = function(id, adj=matrix(0), sign=NULL){
       # extract sign from matrix values if missing sign
       if(is.null(sign)) sign <- as.vector(adj[adj!=0])
       ## create GWASdata object
       new('pathway', id=id, adj=adj, sign=sign)
})

# show method
#' \code{show} displays the \code{\link{pathway}} object briefly   
#' @param object An object of class \code{pathway-class}
## @examples
## #show method
## data(hsa04020)
## hsa04020
#' @export
#' @rdname pathway-class
#' @aliases show,pathway,ANY-method
setMethod('show', signature='pathway',
          definition = function(object){
	      # summarize pathway information
              cat('An object of class ', class(object), ' with id ',object@id,'\n\n',sep='')
              cat('Pathway with adjacency matrix of dimension ', dim(object@adj)[1], ': \n',sep='')
	      # print pathway adjacency matrix
	      net <- object@adj
	      net[net!=0] <- object@sign
              print(net)
              invisible(NULL)
          })

# summary method
setGeneric('summary', function(object, ...) standardGeneric('summary'))

#' \code{summary} generates a \code{\link{pathway}} object summary including basic network properties.
#'
#' @export
#' @rdname pathway-class 
#' @aliases summary,pathway,ANY-method
## @param object An object of class \code{\link{pathway-class}}
## @examples
## #summary method
## data(hsa04020)
## summary(hsa04020)
#' @aliases summary,pathway,ANY-method
setMethod('summary', signature='pathway',
          definition = function(object){
              # define graph and analyze graph
              res <- analyze(object)

              # output
              cat('An object of class ', class(object), '\n\n',sep='')
              cat(res$vcount,' nodes and ',res$ecount,' links; ',res$inh_ecount,' activations and ', res$ecount - res$inh.ecount,' inhibitions. \n\n',sep='')
              cat('Density:\t\t',res$density,' \n')
              cat('Average degree:\t\t',res$av_deg,' \n')
              cat('Inhibition degree:\t',res$inh_deg,' \n')
              cat('Diameter:\t\t',res$diam,' \n')
              cat('Transitivity:\t\t',res$trans,' \n')
              cat('Signed transitivity:\t',res$s_trans,' \n')
              invisible(NULL)
          })

setGeneric('pathway2igraph', function(object, ...) standardGeneric('pathway2igraph'))

#' \code{pathway2igraph} converts a \code{\link{pathway}} object into an 
#' \code{\link[igraph]{igraph}} object with edge attribute \code{sign}
#'
##\link[pkg2:foo_Rd_file_name]{foo}
## @param object An object of class \code{\link{pathway-class}}
#'
#' @return \code{pathway2igraph} returns an unweighted \code{\link[igraph]{igraph}} object with edge attribute \code{sign}
#'
#' @examples
#' # convert to igraph object
#' data(hsa04020)
#' str(hsa04020)
#' g <- pathway2igraph(hsa04020)
#' str(g)
#'
#' @export
#' @rdname pathway-class
#' @aliases pathway2igraph pathway ANY-method     
setMethod('pathway2igraph', signature='pathway', definition = function(object){
    # define adjacency matrix
    net <- object@adj
    net[net!=0] <- object@sign

    # define igraph object
    g <- igraph::graph_from_adjacency_matrix(net, mode='undirected', 
                                     weighted=TRUE, diag=FALSE)

    if(igraph::ecount(g)>0){
      # specify interaction type 
      igraph::E(g)$sign <- igraph::E(g)$weight 
      igraph::E(g)$weight <- 1
    }
    return(g)
})

##### analyze method - analyze pathway network properties
setGeneric('analyze', function(object, ...) standardGeneric('analyze'))

#' analyze \code{\link{pathway}} network properties
#'
#' @export
#' @describeIn pathway
#' 
#' @aliases analyze pathway ANY-method
## @param object An object of class \code{\link{pathway-class}}
#' @return \code{analyze} returns a \code{data.frame} consisting of 
#'   \describe{
#'    \item{id}{pathway id,} 
#'    \item{vcount}{number of genes,}
#'    \item{ecount}{number of links,}
#'    \item{inh_ecount}{number of inhibition links,}
#'    \item{density}{network density,}
#'    \item{av_deg}{average degree,}
#'    \item{inh_deg}{average degree of inhibition links,}
#'    \item{diam}{network diamter,}
#'    \item{trans}{transitivity, and }
#'    \item{s_trans}{signed transitivity (Kunegis et al., 2009).}
#' }
#' @references 
#' Details to the computation and interpretation can be found in:
#' \itemize{
#'   \item Kolaczyk, E. D. (2009). Statistical analysis of network data: methods and models. Springer series in statistics. Springer.
#'   \item Kunegis, J., A. Lommatzsch, and C. Bauckhage (2009). The slashdot zoo: Mining a social network with negative egdes. In Proceedings of the 18th international conference on World wide web, pp. 741-750. ACM Press.
#' } 
#' @examples
#' # analyze pathway network properties
#' data(hsa04020)
#' summary(hsa04020)
#' analyze(hsa04020)
#'
setMethod('analyze', signature='pathway',
          definition = function(object, ...){
              # define graph
              net <- object@adj
              net[net!=0] <- object@sign
              g <- pathway2igraph(object)

              # compute graph characteristics
              res <-  data.frame(id=object@id,
                   vcount=igraph::vcount(g),
                   ecount=igraph::ecount(g),
                   # number of inhibition links
                   inh_ecount = sum(igraph::E(g)$sign<0),
                   # graph density
                   density=igraph::graph.density(g),
                   # average degree
                   av_deg=mean(igraph::degree(g)),
                   # inhibition degree (av neg links)
                   inh_deg=mean(rowSums(net<0)),
                   # diameter
                   diam=igraph::diameter(g),
                   # transitivity
                   trans=igraph::transitivity(g, type='global'),
                   # signed trasitivity (Kunegis, 2009)
                   s_trans = sum(net*(net%*%t(net)))/sum(abs(net)%*%t(abs(net)))
                   )

              # output 
              return(res)
          })

#### function extracting genes in pathway
## @exportMethod get_genes
setGeneric('get_genes', function(object, ...) standardGeneric('get_genes'))

#' \code{get_genes} is a helper function that extracts the gene names in a 
#' \code{\link{pathway}} and returns a \code{vector} containing \code{character}
#' elements of gene names
#'
#' @return \code{get_genes} returns a character vector of gene names extracted from adjacency matrix rownames.
#'
#' @export
#' @describeIn pathway
#' 
#' @aliases get_genes pathway ANY-method 
#' @examples
#' # extract gene names from pathway object
#' get_genes(hsa04020)
#'
setMethod('get_genes', signature='pathway',
          definition = function(object){
              return(rownames(object@adj))
          })

#' @exportMethod plot
if (!isGeneric("plot")) setGeneric('plot')

#' \code{plot} plots \code{\link{pathway}} as \code{\link[igraph]{igraph}} object
#'
#' @export 
#' @rdname pathway-class
#' @aliases plot,pathway,ANY-method
#' 
#' @param x \code{\link{pathway}} object
#' @param y missing (placeholder)
#' @param highlight.genes vector of gene names or node id's, which should be highlighted in a different color, default is \code{NULL} so that no genes are highlighted
#' @param gene.names character indicating whether the genes names should appear in a legend (\code{'legend'}), as vertex label (\code{'nodes'}), or should be omitted (\code{NA})
#' @param main optional overall main title, default is \code{NULL}, which uses the \code{\link{pathway}} id
#' @param asp a \code{numeric} constant, which gives the aspect ratio parameter for plot, default is 0.95
#' @param vertex.size a \code{numeric} constant specifying the vertex size, default is 11
#' @param vertex.color a \code{character} or \code{numeric} constant specifying the vertex color, default is 'khaki1'
#' @param vertex.label.cex a \code{numeric} constant specifying the the vertex label size, default is 0.8,
#' @param edge.width a \code{numeric} constant specifying the edge width, default is 2
#' @param edge.color a \code{character} or \code{numeric} constant specifying the edge color, default is 'olivedrab4'
#' @param ... further arguments specifying plotting options in \code{plot.igraph}
#'
#' @examples
#' # plot pathway as igraph object
#' plot(hsa04020)
#' sample3 <- sample_genes(hsa04020, no = 3)
#' plot(hsa04020, highlight.genes = sample3)
#'
#' @importFrom graphics plot
setMethod('plot', signature(x='pathway',y='missing'),
          function(x, y=NA,
                   highlight.genes = NULL, 
                   gene.names=c('legend','nodes',NA), main=NULL, 
                   asp=0.95, vertex.size=11, vertex.color='khaki1', 
                   vertex.label.cex=0.8,
                   edge.width=2, edge.color='olivedrab4', ...){
                         
    # define igraph object
    g <- pathway2igraph(x)
    
    # define vertex label
    gene.names <- match.arg(gene.names)
    if(gene.names == 'legend')  vertex.label <- 1:length(igraph::V(g))
    if(is.null(gene.names))    vertex.label <- NA
    if(gene.names == 'nodes'){
       vertex.label <- igraph::V(g)$names
    }
    
    # define main title
    if(is.null(main))  main <- paste(x@id)
    
    # interaction type
    igraph::E(g)$lty <- ifelse( igraph::E(g)$sign > 0, 1, 2)
    signs <- ifelse( igraph::E(g)$sign > 0, '+', '-')
    
    # highlighting genes
    if(!is.null(highlight.genes)){
        if(!is.vector(highlight.genes)) stop('highlight.genes has to be a vector')
        vertex.color <- rep(vertex.color, igraph::vcount(g))
        vertex.color[highlight.genes] <- 'yellowgreen' #'darkgreen' 
    }
    
    # define layout
    ltr <- try(igraph::layout.reingold.tilford(g, circular=FALSE), silent=TRUE)
    
    # plot igraph object
    igraph::plot.igraph(g, layout=ltr, asp=asp, 
                vertex.size=vertex.size, vertex.color=vertex.color, vertex.label.cex=vertex.label.cex, vertex.label=vertex.label, 
                edge.width=edge.width, edge.color=edge.color,
                main=main, ...)
    # add legend
    if(gene.names == 'legend'){
      graphics::legend(,x=-1.55, y=1.5, bty='n',lty=NULL, legend=paste(1:igraph::vcount(g),igraph::V(g)$name),cex=0.7)
    }
    
    invisible(NULL)
})

#' @exportMethod sample_genes
setGeneric('sample_genes', function(object, ...) standardGeneric('sample_genes'))

#' \code{sample_genes} function randomly selects effect gene in a 
#' \code{\link{pathway}} according the betweenness centrality and (no -1) neighors
#'
#' @return \code{sample_genes} returns a \code{vector} of length \code{no} with 
#' vertex id's of sampled genes 
#'
#' @export
#' @describeIn pathway
#'
#' @param no a \code{numeric} constant specifying the number of genes to be sampled, default is 3
#' @aliases sample_genes pathway ANY-method 
#' @examples
#' # sample effect genes
#' sample3 <- sample_genes(hsa04020, no = 3)
#' plot(hsa04020, highlight.genes = sample3)
#' sample5 <- sample_genes(hsa04020, no = 5)
#' plot(hsa04020, highlight.genes = sample5)
#'
setMethod('sample_genes', signature='pathway',
          definition = function(object, no=3){
            g <- pathway2igraph(object)
            if(igraph::vcount(g) < no) stop('number of sampled genes should be smaller 
                               than total number of genes in the pathway')

            # sample a gene with high betweennes centrality
            sel1 <- sample(igraph::V(g), size=1, prob=igraph::betweenness(g))
            # sample two of its neighbors
            sel2 <- sample(igraph::ego(g,1,sel1)[[1]], size=no-1)
            # combine samples
            samp <- c(as.numeric(sel1),sel2)
            names(samp)[1] <- sel1$name

            return(samp)
          })


setGeneric('gene_name_number', function(x, ...) standardGeneric('gene_name_number'))
#' Function to get genes names and numbers from kegg (for internal use)
#'
#' This function extracts for a particular \code{\link{pathway}} all genes and 
#' the numbers they are represented with in the KEGG network from the 
#' corresponding KGML pathway file.
#'
#' @param x A \code{character} hsa identifier of the pathway for which gene 
#' infomation should be extracted as used in KEGG database ('hsa00100').    
#' @return A \code{data.frame} listing the genes included in the pathway with 
#' their names as well as numbers used in KEGG database.
#' @author Stefanie Friedrichs, Patricia Burger
#' @keywords internal
setMethod('gene_name_number', signature='character', 
          definition = function(x){
    info <- scan(url(paste("http://togows.dbcls.jp/entry/pathway/",
                           x,"/genes",sep="")), what="character")   
    pos <- which(substr(info,nchar(info),nchar(info))==";")
    #if ";" entry is pos and pos-2 has "]" as last character -> ok. 
    #else: search until pos-j ends with "]".
    liste     <- matrix(rep(0,length(pos)*2),ncol=2, byrow=TRUE)
    #Above "if else" is not applicable for the first pos. Therefore:
    #If item before first gene name begins not with at least 2 digits search
    #all entries before the first pos for matching pattern
    if(grepl("^[0-9]{2,}" ,info[pos[1]-1]) == FALSE){
      if(pos[1] <= 2){
        cat("First gene in this pathway has incorrect hsa number. Please
            check the KEGG database for this pathway.")
      }      
      for(j in 2:pos[1]-1){
        check <- grepl("^[0-9]{2,}" ,info[pos[1]-j])
        if(check){
          liste[1,] <- info[c((pos[1]-j),pos[1])]
          break
        }
      }
    } else{
      liste[1,] <- info[c((pos[1]-1),pos[1])] #before first entry no "]"
    }
    
    
    if(length(pos) == 1){
      warning("This pathway contains only one gene.")
      print("This pathway contains only one gene.")
      return(liste)
    }
    
    for(i in 2:length(pos)){ 
      if(substr(info[pos[i]-2],nchar(info[pos[i]-2]),nchar(info[pos[i]-2]))=="]"){ 
         liste[i,] <- info[c((pos[i]-1),pos[i])]
         }else{  
           j <- 3 
           while(liste[i,1]=="0"){
           if(substr(info[pos[i]-j],nchar(info[pos[i]-j]),nchar(info[pos[i]-j]))=="]"){
              textt <- paste(info[ (pos[i]-j+2):pos[i] ], collapse = '')
              liste[i,] <- c(info[c(pos[i]-j+1)],textt) 
           }else{ 
             if(grepl(";", info[pos[i]-j])){
               print("This gene has no [KO:XXXXXXX] number.")
               if(grepl("^[0-9]{2,}" ,info[pos[1]-1])){
                 liste[i,] <- info[c((pos[i]-1),pos[i])]
                 
               }
             }
             j <- j+1
           }      
           }
       }
    }
    liste[,2] <- substr(liste[,2],1,nchar(liste[,2])-1) #cut ";"
    return(liste)
})


#' An S4 class for an object assigning genes to pathways
#'
#' @rdname pathway_info-class
#' @slot info A \code{data.frame} including information on genes contained in 
#' pathways with columns 'pathway', 'gene_start', 'gene_end', 'chr' and 'gene'.
#'
#' @author Stefanie Friedrichs, Juliane Manitz
pathway_info <- setClass('pathway_info', slots=c(info='data.frame'))

setValidity('pathway_info', function(object){  
	msg  <- NULL
	valid <- TRUE
 	if(!is.data.frame(object@info)){
	  valid=FALSE
	  msg <- c(msg, "the pathway_info object must include a data.frame")
	}
	if(!all.equal(colnames(object@info),c("pathway","gene_start","gene_end","chr","gene") )){
	  valid=FALSE
	  msg <- c(msg, "the included data.frame needs columns 'pathway',
           'gene_start', 'gene_end', 'chr' and 'gene'.")
	}
})

setGeneric('pathway_info', function(x) standardGeneric('pathway_info'))
#' Get information on genes in a pathway
#'
#' This function lists all genes formig a particular \code{\link{pathway}}. Start and end  
#' positions of these genes are extracted from the Ensemble database. The 
#' database is accessed via the R-package \pkg{biomaRt}.
#'
#' @param x A \code{character} identifying the pathway for which gene infomation 
#' should be extracted. Here KEGG IDs (format: 'hsa00100') are used.
#' @return A \code{data.frame} including as many rows as genes appear in the 
#' \code{\link{pathway}}. for each gene its name, the start and end point and the chromosome 
#' it lies on are given.
#'
#' @import biomaRt  
#' @export
#' @seealso \code{\link{snp_info}}, \code{\link{get_anno}}
#' @rdname pathway_info-class
setMethod('pathway_info', signature='character', 
         definition = function(x){              
   g       <- gene_name_number(x)[,2] 
   #ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")  
   ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
              dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org") 
   info    <- stats::na.omit(getBM(attributes=c("start_position","end_position",
              "chromosome_name","hgnc_symbol"), filters=c("hgnc_symbol"),
              values=g, mart=ensembl))
#   info$chromosome_name <- as.numeric(info$chromosome_name)
#   info                 <- stats::na.omit(info)
#   pathways             <- cbind(rep(paste(x,sep=""),length(info[,1])),info)
   pathways             <- data.frame(pathway=x, info)
   colnames(pathways)   <- c("pathway","gene_start","gene_end","chr","gene") 
   ret <- new('pathway_info', info=pathways)
   return(ret)
})
#' \code{show} Shows basic information on \code{\link{pathway_info}} object
#'
#' @param object An object of class \code{\link{pathway_info}}.
#' @return \code{show} Basic information on \code{\link{pathway_info}} object.
#' @examples
#' data(hsa04022_info) # pathway_info('hsa04020') 
#' show(hsa04022_info)
#' summary(hsa04022_info)
#'
#' @export
#' @rdname pathway_info-class
#' @aliases show,pathway_info,ANY-method
setMethod('show', signature='pathway_info',
          definition = function(object){
            cat('An object of class ', class(object), '\n', sep='')
            cat('Number of pathways:', length(unique(object@info$pathway)),'\n')
            cat('Number of genes:', length(unique(object@info$gene)),'\n')
            if(nrow(object@info)>6){ cat('First six rows: \n')
            print(object@info[1:6,])}else{print(object@info)}   
          })

setGeneric('summary', function(object, ...) standardGeneric('summary'))
#' \code{summary} Summarizes information on \code{\link{pathway_info}} object
#'
## @param object An object of class \code{\link{pathway_info}}.
#' @return \code{summary} Summarized information on \code{\link{pathway_info}} object.
#' @export
#' @rdname pathway_info-class
#' @aliases summary,pathway_info,ANY-method
setMethod('summary', signature='pathway_info',
          definition = function(object){
            cat('An object of class ', class(object), '\n', sep='')
            cat('Number of pathways:', length(unique(object@info$pathway)),'\n')
            cat('Number of genes:', length(unique(object@info$gene)),'\n')
            if(nrow(object@info)>6){ cat('First six rows: \n')
            print(object@info[1:6,])}else{print(object@info)}   
          })


setGeneric('set_one', function(x, ...) standardGeneric('set_one'))
#' Helper function to set \code{matrix} entries to 0/1/-1 only 
#'
#' This function sets all entries in a \code{matrix} bigger than 1 to 1 and all entries 
#' smaller than -1 to -1. It is called by \code{\link{get_network_matrix}}.
#' (For internal use)
#'
#' @param x numeric \code{matrix} representing the network adjacency matrix. 
#' @return A \code{matrix} representing the interaction network in the 
#' \code{\link{pathway}} object with entries equal to 1, -1 or 0.
#'
#' @author Stefanie Friedrichs 
#' @seealso \code{\link{get_network_matrix}}
#' @keywords internal
setMethod('set_one', signature='matrix', 
          definition = function(x){       
  if(length(x[x>1])>0){ 
     print("Edge values greater than 1 set to 1!")
     x[x>1] <- 1}
  if(length(x[x < (-1)])>0){
     print("Edge values smaller than -1 set to -1!")
     x[x<(-1)] <- -1}
  return(x)
})


setGeneric('set_names', function(x, ...) standardGeneric('set_names'))
#' Helper function to translate gene numbers to names  
#'
#' This function exchanges the numbers used for genes in KEGG download KGML files
#' with the corresponding gene names. Names are set to be the column names and 
#' rownames of the pathway's network matrix. The function 
#' is called by \code{\link{get_network_matrix}}. (For internal use)
#'
#' @param x A \code{matrix} representing the network matrix. 
#' @param nodes A \code{vector} of gene numbers to be replaced by names.   
#' @param my_list A \code{data.frame} listing gene names and numbers. Output from 
#' \code{gene_name_number}.  
#' @return A \code{matrix} representing the interaction network in the pathway
#' with gene names as rownames and columnnames. 
#'
#' @author Stefanie Friedrichs
#' @keywords internal
setMethod('set_names', signature='matrix', 
          definition = function(x, nodes, my_list){
    name <- substr(nodes,5,nchar(nodes)) 
    for(i in 1:length(name)){ 
      name[i] <- my_list[my_list[,1]==name[i],2] 
    }
    colnames(x) <- rownames(x) <- name
    return(x)
})
 
 
setGeneric('get_network_matrix', function(object, ...) standardGeneric('get_network_matrix'))
#' Function to calculate the network \code{matrix} for a \code{\link{pathway}} object
#'
#' \code{get_network_matrix} creates the adjacency matrix representing the gene-gene interaction structure within a particular \code{\link{pathway}}. Note that a 
#' KEGG kgml file is downloaded and saved in the working directory. 
#'
#' @param object A \code{\link{pathway}} object identifying the pathway for which gene interaction infomation should be extracted. Here, KEGG IDs of format 'hsa00100' are used and information is downloaded from the KEGG database.    
#' @param directed A \code{logical} argument, stating whether the network matrix 
#' should return directed (\code{TRUE}) or undirected (\code{FALSE}) links. 
#' @param method Download method to be used for downloading files, passed to via \code{KEGGgraph::retrieveKGML} to \code{utils::download.file} function. Currently supports \code{'auto'} (default), \code{'internal'}, \code{'wininet'} (Windows only), \code{'libcurl'}, \code{'wget'} and \code{'curl'}. 
#' @return \code{get_network_matrix} returns the modified \code{\link{pathway}} object, where the slots \code{adj} and \code{sign} are altered according to the downloaded information in the KEGG kgml file. 
#' @examples
#' get_network_matrix(pathway(id="hsa04020"), directed=TRUE)
#'
#' @author Stefanie Friedrichs, Patricia Burger, Juliane Manitz
#' @export 
#' @name get_network_matrix
#' @rdname get_network_matrix-methods
#' @aliases get_network_matrix,pathway-method
#' @importFrom KEGGgraph retrieveKGML
#' @importFrom KEGGgraph parseKGML2Graph
#' @importFrom KEGGgraph nodes
#' @importFrom KEGGgraph parseKGML2DataFrame
setMethod('get_network_matrix', signature='pathway', 
          definition = function(object, directed=TRUE, method='auto'){    
      
      x <- object@id      
      retrieveKGML(substr(x,4,nchar(x)), organism="hsa",
                   destfile=paste(x,".xml",sep=""), method=method)
      liste     <- gene_name_number(x)
      pathgraph <- parseKGML2Graph(paste(x,".xml",sep=""), expandGenes=TRUE)
          
     # parseKGML2Graph can also be split up in two steps 
     # pathgraph.s1 <- parseKGML(paste(x, ".xml", sep = ""))
     # pathgraph  <- KEGGpathway2Graph(pathgraph.s1, expandGenes = TRUE)
           
      nodes     <- nodes(pathgraph) #vector of gene no. (format: "hsa:226")
      edgelist  <- parseKGML2DataFrame(paste(x,".xml",sep=""))
      if(length(edgelist) != 0){
        edgelist  <- edgelist[!is.na(edgelist[,3]),] #delete NA-types 
      }
         
      # --- Skip pathway, if no edges or wrong number of genes ---
      if(length(edgelist)==0){
        print(paste("KGML-file for ",x," has no edges!",sep=""))
        }   
      if(length(nodes)!=length(liste[,2])){
        stop(paste("Wrong number of genes in ",x,"!",sep=""))
        }               
      #set up empty matrix:
      N <- matrix(rep(0, (length(nodes)^2)),nrow=length(nodes))
      colnames(N) <- rownames(N) <- nodes 
      verb.a <- NULL
      verb.i <- NULL
      verb.s <- NULL
      if(length(edgelist) != 0){
        verb.a <- edgelist[edgelist[,3]=="activation",] 
        verb.i <- edgelist[edgelist[,3]=="inhibition",]
        verb.s <- edgelist[edgelist[,3]!="inhibition"&edgelist[,3]!="activation",] 
      }
             
      # --- if activations/inhibitions specified: Signed matrix ---
      if( (length(verb.a[,1])+length(verb.i[,1]))>0 ){
        if(length(verb.s[,1])==0){ 
          print(paste(x,": Activation/Inhibition only: 
                      Signed graph!",sep="")) 
        }
        if(length(verb.s[,1])>0){ 
          print(paste(x," has both: Activation/Inhibition edges and edges without type!",sep=""))  
                      
        }
        # -- Directed --- 
        #           to 
        # from  (        )
        if(length(verb.a[,1])>0){ 
          from <- as.character(unique(verb.a[,1]))
          for(i in from){ 
            N[i,as.character(verb.a[verb.a[,1]==i,2])] <- 1}  
        }
        if(length(verb.i[,1])>0){ 
          from <- as.character(unique(verb.i[,1]))
          for(i in from){ 
            if(sum(N[i,as.character(verb.i[verb.i[,1]==i,2])])>0){
              print("Edge will be removed!")} #if i. at same edge as a. was
            N[i,as.character(verb.i[verb.i[,1]==i,2])] <- 
              N[i,as.character(verb.i[verb.i[,1]==i,2])] -1 }   } 
        N <- set_names(N,nodes,liste)                    
        # -- Undirected --- 
        M <- N + t(N) #contradictory edges will be removed
        M <- set_one(M)    #double edges (2,-2) could be produced
        M <- set_names(M,nodes,liste)           
        }
      
      # --- else, if no edge types specified: Unsigned matrix ---  
      #only use edges without type, if ALL edges are without type
      if(length(verb.s[,1])>0 & (length(verb.a[,1])+length(verb.i[,1]))==0){
        print(paste(x,": No edge-types, unsigned graph only!",sep=""))
        
        # -- Directed ---
        if(length(verb.s[,1])>0){ 
          from <- as.character(unique(verb.s[,1]))
          for(i in from){ 
            N[i,as.character(verb.s[verb.s[,1]==i,2])] <- 1
          }  
        }
        N <- set_names(N,nodes,liste)                           
        # -- Undirected --- 
        M <- N + t(N)
        #contradictory edges will be removed
        #double edges (2,-2) can be produced:
        M <- set_one(M)
        M <- set_names(M,nodes,liste)      
      }
      if(directed==TRUE){
         object@adj  <- abs(N) 
         object@sign <- as.vector(N[N!=0])
         return(object)
      }else{
         object@adj  <- abs(M) 
         object@sign <- as.vector(M[M!=0])
         return(object)
      }                       
})
