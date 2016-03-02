#################################################
#
# pathway object functions
#
#################################################

#' An S4 class to represent a gene-gene interaction network
#'
#' @rdname pathway-class
#'
#' @slot id a character repesenting the pathway id, e.g. hsa00100 as used in 
#' the KEGG database.
#' @slot adj a matrix respresenting the network adjacency matrix  of dimension 
#' equaling the number of genes (1 interaction, 0 otherwise)
#' @slot sign a numeric vector indicating the interaction type for each link 
#' (1 activation, -1 inhibition). Represents the interaction network within the 
#' pathway.
#' 
#' @author Juliane Manitz
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

# show method
#' \code{show} displays the pathway object briefly
#' @param object pathway object
#' @examples
#' #show method
#' data(hsa04020)
#' hsa04020
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

#' \code{summary} generates a pathway object summary including basic network properties.
#'
#' @export
#' @rdname pathway-class 
#' @aliases summary,pathway,ANY-method
#'
#' @examples
#' #summary method
#' data(hsa04020)
#' summary(hsa04020)
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

#' \code{pathway2igraph} converts a \code{\link{pathway}} object into an \code{\link{igraph}} object with edge attribute \code{sign}
#'
#' @export
#' @rdname pathway-class 
#' @aliases pathway2igraph,pathway,ANY-method
#' 
#' @examples
#' # convert to igraph object
#' data(hsa04020)
#' str(hsa04020)
#' g <- pathway2igraph(hsa04020)
#' str(g)
#'
#' @import igraph
setMethod('pathway2igraph', signature='pathway',
          definition = function(object){
    # define adjacency matrix
    net <- object@adj
    net[net!=0] <- object@sign

    # define igraph object

    g <- graph_from_adjacency_matrix(net, mode='undirected', 
                                     weighted=TRUE, diag=FALSE)

    if(ecount(g)>0){
      # specify interaction type 
      E(g)$sign <- E(g)$weight 
      E(g)$weight <- 1
    }
    return(g)
})

##### analyze method - analyze pathway network properties
setGeneric('analyze', function(object, ...) standardGeneric('analyze'))

#' analyze pathway network properties
#'
#' @export
#' @rdname pathway-class
#' @aliases analyze,pathway,ANY-method
#'
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
#' data(hsa04020)
#' summary(hsa04020)
#' analyze(hsa04020)
#'
#' @import igraph
setMethod('analyze', signature='pathway',
          definition = function(object, ...){
              # define graph
              net <- object@adj
              net[net!=0] <- object@sign
              g <- pathway2igraph(object)

              # compute graph characteristics
              res <-  data.frame(id=object@id,
                   vcount=vcount(g),
                   ecount=ecount(g),
                   # number of inhibition links
                   inh_ecount = sum(E(g)$sign<0),
                   # graph density
                   density=graph.density(g),
                   # average degree
                   av_deg=mean(degree(g)),
                   # inhibition degree (av neg links)
                   inh_deg=mean(rowSums(net<0)),
                   # diameter
                   diam=diameter(g),
                   # transitivity
                   trans=transitivity(g, type='global'),
                   # signed trasitivity (Kunegis, 2009)
                   s_trans = sum(net*(net%*%t(net)))/sum(abs(net)%*%t(abs(net)))
                   )

              # output 
              return(res)
          })

#### function extracting genes in pathway
#' @exportMethod get_genes
setGeneric('get_genes', function(object, ...) standardGeneric('get_genes'))

#' \code{get_genes} is a helper function that extracts the gene names in a pathway and returns a vector of character containing gene names
#'
#' @export
#' @rdname pathway-class 
#' @aliases get_genes,pathway,ANY-method
#' @param object A pathway object
#' @examples
#' # extract gene names from pathway
#' get_genes(hsa04020)
setMethod('get_genes', signature='pathway',
          definition = function(object){
              return(rownames(object@adj))
          })

#' @exportMethod plot
if (!isGeneric("plot")) setGeneric('plot')

#' \code{plot} plots pathway as igraph object
#'
#' @export 
#' @rdname pathway-class
#' @aliases plot,pathway,ANY-method
#' 
#' @param x pathway object
#' @param y missing (placeholder)
#' @param highlight.genes vector of gene names or node id's, which should be highlighted in a different color, default is \code{NULL} so that no genes are highlighted
#' @param gene.names character indicating whether the genes names should appear in a legend (\code{'legend'}), as vertex label (\code{'nodes'}), or should be omitted (\code{NA})
#' @param main optional overall main title, default is \code{NULL}, which uses the pathway id
#' @param asp a numeric constant, which gives the aspect ratio parameter for plot, default is 0.95
#' @param vertex.size a numeric constant specifying the vertex size, default is 11
#' @param vertex.color a character or numeric constant specifying the vertex color, default is 'khaki1'
#' @param vertex.label.cex a numeric constant specifying the the vertex label size, default is 0.8,
#' @param edge.width a numeric constant specifying the edge width, default is 2
#' @param edge.color a character or numeric constant specifying the edge color, default is 'olivedrab4'
#' @param ... further arguments specifying plotting options in \code{\link{plot.igraph}}
#'
#' @examples
#' # plot pathway as igraph object
#' plot(hsa04020)
#' sample3 <- sample_genes(hsa04020, no = 3)
#' plot(hsa04020, highlight.genes = sample3)
#'
#' @import igraph
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
if(gene.names == 'legend')  vertex.label <- 1:length(V(g))
if(is.null(gene.names))    vertex.label <- NA
if(gene.names == 'nodes'){
   vertex.label <- V(g)$names
}

# define main title
if(is.null(main))  main <- paste(x@id)

# interaction type
E(g)$lty <- ifelse( E(g)$sign > 0, 1, 2)
signs <- ifelse( E(g)$sign > 0, '+', '-')

# highlighting genes
if(!is.null(highlight.genes)){
    if(!is.vector(highlight.genes)) stop('highlight.genes has to be a vector')
    vertex.color <- rep(vertex.color, vcount(g))
    vertex.color[highlight.genes] <- 'yellowgreen' #'darkgreen' 
}

# define layout
ltr <- try(layout.reingold.tilford(g, circular=FALSE), silent=TRUE)

# plot igraph object
plot.igraph(g, layout=ltr, asp=asp, 
            vertex.size=vertex.size, vertex.color=vertex.color, vertex.label.cex=vertex.label.cex, vertex.label=vertex.label, 
            edge.width=edge.width, edge.color=edge.color,
            main=main, ...)
# add legend
if(gene.names == 'legend'){
  legend(,x=-1.55, y=1.5, bty='n',lty=NULL, legend=paste(1:vcount(g),V(g)$name),cex=0.7)
}

invisible(NULL)
})

#' @exportMethod sample_genes
setGeneric('sample_genes', function(object, ...) standardGeneric('sample_genes'))

#' \code{sample_genes} function randomly selects effect genes in pathway and returns a vector of length \code{no} with vertex id's of sampled genes 
#'
#' @export
#' @rdname pathway-class 
#' @aliases get_genes,pathway,ANY-method
#'
#' @param no a numeric constant specifying the number of genes to be sampled, default is 3
#' 
#' @examples
#' # sample effect genes
#' sample3 <- sample_genes(hsa04020, no = 3)
#' plot(hsa04020, highlight.genes = sample3)
#' sample5 <- sample_genes(hsa04020, no = 5)
#' plot(hsa04020, highlight.genes = sample5)
#'
#' @import igraph
setMethod('sample_genes', signature='pathway',
          definition = function(object, no=3){
            g <- pathway2igraph(object)
            if(vcount(g) < no) stop('number of sampled genes should be smaller 
                               than the total number in the pathway')

            # sample a gene with high betweennes centrality
            sel1 <- sample(V(g), size=1, prob=betweenness(g))
            # sample two of its neighbors
            sel2 <- sample(ego(g,1,sel1)[[1]], size=no-1)
            # combine samples
            samp <- c(as.numeric(sel1),sel2)
            names(samp)[1] <- sel1$name

            return(samp)
          })


setGeneric('gene_name_number', function(x, ...) standardGeneric('gene_name_number'))
#' Function to get genes names and numbers from kegg (for internal use)
#'
#' This function extracts for a particular pathway all genes and the numbers
#' they are represented with in the KEGG network from the corresponding KGML pathway file.
#'
#' @param x A \code{character} hsa identifier of the pathway for which gene 
#' infomation should be extracted as used in KEGG database ('hsa00100').    
#' @return A \code{data.frame} listing the genes included in the pathway with 
#' their names as well as numbers used in KEGG database.
#' @author Stefanie Friedrichs
#' @keywords internal
setMethod('gene_name_number', signature='character', 
          definition = function(x){
    info <- scan(url(paste("http://togows.dbcls.jp/entry/pathway/",
                           x,"/genes",sep="")), what="character")   
    pos <- which(substr(info,nchar(info),nchar(info))==";")
    #if ";" entry is pos and pos-2 has "]" as last character -> ok. 
    #else: search until pos-j ends with "]".
    liste     <- matrix(rep(0,length(pos)*2),ncol=2, byrow=TRUE)
    liste[1,] <- info[c((pos[1]-1),pos[1])] #before first entry no "]"
    for(i in 2:length(pos)){ 
      if(substr(info[pos[i]-2],nchar(info[pos[i]-2]),nchar(info[pos[i]-2]))=="]"){ 
         liste[i,] <- info[c((pos[i]-1),pos[i])]}else{  
           j <- 3 
           while(liste[i,1]=="0"){
           if(substr(info[pos[i]-j],nchar(info[pos[i]-j]),nchar(info[pos[i]-j]))=="]"){
              textt <- paste(info[ (pos[i]-j+2):pos[i] ], collapse = '')
              liste[i,] <- c(info[c(pos[i]-j+1)],textt) }else{ j <- j+1 }      
           }
       }
    }
    liste[,2] <- substr(liste[,2],1,nchar(liste[,2])-1) #cut ";"
    return(liste)
})


#' An S4 class for an object assigning genes to pathways
#'
#' @rdname spathway_info-class
#' @slot info A \code{data.frame} including information on genes contained in 
#' pathways
#' @examples
#' data(hsa04022_info) 
#'
#' @author Stefanie Friedrichs
#' @export pathway_info
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

setGeneric('pathway_info', function(x, ...) standardGeneric('pathway_info'))
#' Get information on genes in a pathway
#'
#' This function lists all genes formig a particular pathway. Start and end  
#' positions of these genes are extracted from the Ensemble database. The 
#' database is accessed via the R-package \code{biomaRt}.
#'
#' @param x A \code{character} identifying the pathway for which gene infomation 
#' should be extracted. Here KEGG IDs ('hsa00100') are used.
#' @return A \code{data.frame} including as many rows as genes appear in the 
#' pathway. for each gene its name, the start and end point and the chromosome 
#' it lies on are given.
#' @examples
#' pathway_info("hsa04022") 
#'
#' @author Stefanie Friedrichs
#' @import biomaRt  
#' @export
#' @rdname pathway_info-class
setMethod('pathway_info', signature='character', 
         definition = function(x){              
   g       <- gene_name_number(x)[,2] 
   #ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")  
   ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
              dataset="hsapiens_gene_ensembl",host = "jul2015.archive.ensembl.org") 
   info    <- getBM(attributes=c("start_position","end_position",
              "chromosome_name","hgnc_symbol"), filters=c("hgnc_symbol"),
              values=g, mart=ensembl)
   info$chromosome_name <- as.numeric(as.character(info$chromosome_name))
   info                 <- na.omit(info)
   pathways             <- cbind(rep(paste(x,sep=""),length(info[,1])),info)
   colnames(pathways)   <- c("pathway","gene_start","gene_end","chr","gene") 
   ret <- new('pathway_info', info=pathways)
   return(ret)
})
#' \code{show} Shows basic information on \code{pathway_info} object
#'
#' @param object An object of class \code{\link{pathway_info}}.
#' @return \code{show} Basic information on \code{pathway_info} object.
## @author Stefanie Friedrichs
#' @examples
#' # show method
#' data(hsa04022_info)
#' hsa04022_info
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
#' \code{summary} Summarizes information on \code{pathway_info} object
#'
#' @param object An object of class \code{\link{pathway_info}}.
#' @return \code{summary} Summarized information on \code{pathway_info} object.
## @author Stefanie Friedrichs
#' @examples
#' # summary method
#' data(hsa04022_info)
#' summary(hsa04022_info)
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
#' Helper function to set matrix entries to 0/1/-1 only 
#'
#' This function sets all entries in a matrix bigger than 1 to 1 and all entries 
#' smaller than -1 to -1. It is called by \code{\link{get_network_matrix}}.
#' (For internal use)
#'
#' @param x numeric \code{matrix} representing the network adjacency matrix. 
#' @return A \code{matrix} representing the interaction network in the pathway
#' with entries equal to 1, -1 or 0.
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
#' @param x A \code{matrix} representing the networkmatrix. 
#' @param nodes A \code{vector} of gene numbers to be replaced by names.   
#' @param my_list A \code{data.frame} listing gene names and numbers. Output from 
#' \code{gene_name_number}.  
#' @return A \code{matrix} representing the interaction network in the pathway
#' with gene names as rownames and columnnames. 
#'
#' @author Stefanie Friedrichs
#' @seealso \code{\link{get_network_matrix}}
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
 
 
setGeneric('get_network_matrix', function(x, ...) standardGeneric('get_network_matrix'))
#' Function to calculate the network matrix
#'
#' This function creates the networkmatrix representing the gene-gene interaction 
#' structure within a particular pathway. In this process a KEGG kgml file is 
#' downloaded and saved in the working directory. 
#'
#' @param x A \code{character} identifying the pathway for which gene infomation 
#' should be extracted. Here KEGG IDs ('hsa00100') are used. 
#' @param directed A \code{logic} argument, stating whether the networkmatrix 
#' should be returned directed (\code{TRUE}) or undirected (\code{FALSE}).
#' @param keep.kgml A \code{logic} argument, specifying whether the downloaded 
#' KEGG kgml file of the pathway should be kept in the working directory. 
#' For (\code{FALSE}) the file is deleted, for (\code{TRUE}) not.  
#' after calculation of the network matrix. 
#' @return A \code{matrix} representing the interaction network in the pathway.
## @examples
## get_network_matrix("hsa04022", TRUE)
#'
#' @author Stefanie Friedrichs
#' @import KEGGgraph, biomaRt 
#' @export   
setMethod('get_network_matrix', signature='character', 
          definition = function(x, directed, keep.kgml){    
    retrieveKGML(substr(x,4,nchar(x)), organism="hsa",
                 destfile=paste(x,".xml",sep=""), method="internal")
    #filename  <- paste(x,".xml",sep="")
    liste     <- gene_name_number(x)
    pathgraph <- parseKGML2Graph(paste(x,".xml",sep=""), expandGenes=TRUE)
    nodes     <- nodes(pathgraph) #vector with gene numbers (format: "hsa:226")
    edgelist  <- parseKGML2DataFrame(paste(x,".xml",sep=""))
    edgelist  <- edgelist[!is.na(edgelist[,3]),] #delete NA-types            

  # --- Skip pathway, if no edges or wrong number of genes ---
    if(length(edgelist)==0){
       print(paste("KGML-file for ",x," has no edges!",sep=""))
       next}   
    if(length(nodes)!=length(liste[,2])){
       paste(paste("Wrong number of genes in ",x,"!",sep=""))
       next}               
    #set up empty matrix:
    N <- matrix(rep(0, (length(nodes)^2)),nrow=length(nodes))
    colnames(N) <- rownames(N) <- nodes 
    verb.a <- edgelist[edgelist[,3]=="activation",] 
    verb.i <- edgelist[edgelist[,3]=="inhibition",]
    verb.s <- edgelist[edgelist[,3]!="inhibition"&edgelist[,3]!="activation",] 
  
  # --- if activations/inhibitions specified: Signed matrix ---
    if( (length(verb.a[,1])+length(verb.i[,1]))>0 ){
      if(length(verb.s[,1])==0){ print(paste(x,": Activation/Inhibition only: 
                                 Signed graph!",sep="")) }
      if(length(verb.s[,1])>0){ print(paste(x," has both: Activation/Inhibition   
                                      edges and edges without type!",sep=""))}
     # -- Directed --- 
     #           to 
     # from  (        )
     if(length(verb.a[,1])>0){ 
        from <- as.character(unique(verb.a[,1]))
        for(i in from){ 
           N[i,as.character(verb.a[verb.a[,1]==i,2])] <- 1}  }
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
    if(keep.kgml == FALSE){ 
       rm(paste(x,".xml",sep="")) 
    }
    if(directed==TRUE){
      return(N)
    }else{
      return(M)
    }    
})
          
