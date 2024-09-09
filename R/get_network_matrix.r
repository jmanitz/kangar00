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
# #' @examples
# #' get_network_matrix(pathway(id="hsa04020"), directed=TRUE)
#'
#' @author Stefanie Friedrichs, Patricia Burger, Juliane Manitz
#' @export 
#' @name get_network_matrix
#' @rdname get_network_matrix
#' @aliases get_network_matrix,pathway-method
# @importFrom KEGGgraph retrieveKGML
# @importFrom KEGGgraph parseKGML2Graph
# @importFrom KEGGgraph nodes
# @importFrom KEGGgraph parseKGML2DataFrame
setMethod('get_network_matrix', signature='pathway', 
          definition = function(object, directed=TRUE, method='auto'){    
            
            x <- object@id      
            KEGGgraph::retrieveKGML(substr(x,4,nchar(x)), organism="hsa",
                         destfile=paste(x,".xml",sep=""), method=method)
            liste     <- gene_name_number(x)
            pathgraph <- KEGGgraph::parseKGML2Graph(paste(x,".xml",sep=""), expandGenes=TRUE)
            
            # parseKGML2Graph can also be split up in two steps 
            # pathgraph.s1 <- parseKGML(paste(x, ".xml", sep = ""))
            # pathgraph  <- KEGGpathway2Graph(pathgraph.s1, expandGenes = TRUE)
            
            nodes     <- KEGGgraph::nodes(pathgraph) #vector of gene no. (format: "hsa:226")
            edgelist  <- KEGGgraph::parseKGML2DataFrame(paste(x,".xml",sep=""))
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

