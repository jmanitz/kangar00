#################################################
#
# pathway object functions
#
#################################################

# object constructor 
pathway <- setClass('pathway',
                    slots=c(id='character', adj='matrix', sign='vector'))

# id .. pathway id, e.g. hsa00100
# adj .. adjacency matrix
# sign .. interaction type for each link

# validy checks
setValidity('pathway', function(object){ 
	msg  <- NULL
	valid <- TRUE
	# adjacency includes only 1 and 0
	if(!all( object@adj %in% c(0,1))){
	  valid=FALSE
	  msg <- c(msg, "the adjacency matrix is not allowed to include other values that zero and one")
	}
	# length(sign) = number of links
	if(length(object@sign)!=sum(object@adj!=0)){
	  valid <- FALSE
          msg <- c(msg, "the length of the sign vector does not correspond to the number of links")
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

# convert pathway to igraph object
pathway2igraph <- function(pathway){

    # define adjacency matrix
    net <- pathway@adj
    net[net!=0] <- pathway@sign

    # define igraph object
    g <- graph.adjacency(net, mode='undirected', weighted=TRUE, diag=FALSE)

    # specify interaction type 
    E(g)$sign <- E(g)$weight 
    E(g)$weight <- 1

    return(g)
}

# show method
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

setMethod('summary', signature='pathway',
          definition = function(object){
              # define graph and analyze graph
              res <- analyze(object)

              # output
              cat('An object of class ', class(object), '\n\n',sep='')
              cat(res$vcount,' nodes and ',res$ecount,' links; ',res$inh.ecount,' activations and ', res$ecount - res$inh.ecount,' inhibitions. \n\n',sep='')
              cat('Density:\t\t',res$density,' \n')
              cat('Average degree:\t\t',res$av.deg,' \n')
              cat('Inhibition degree:\t',res$inh.deg,' \n')
              cat('Diameter:\t\t',res$diam,' \n')
              cat('Transitivity:\t\t',res$trans,' \n')
              cat('Signed transitivity:\t',res$s.trans,' \n')
              invisible(NULL)
          })

# analyze method - analyze pathway network properties
setGeneric('analyze', function(object, ...) standardGeneric('analyze'))

setMethod('analyze', signature='pathway',
          definition = function(object){
              # define graph
              net <- object@adj
              net[net!=0] <- object@sign
              g <- pathway2igraph(object)

              # compute graph characteristics
              res <-  data.frame(id=object@id,
                   vcount=vcount(g),
                   ecount=ecount(g),
                   # number of inhibition links
                   inh.ecount = sum(E(g)$sign<0),
                   # graph density
                   density=graph.density(g),
                   # average degree
                   av.deg=mean(degree(g)),
                   # inhibition degree (av neg links)
                   inh.deg=mean(rowSums(net<0)),
                   # diameter
                   diam=diameter(g),
                   # transitivity
                   trans=transitivity(g, type='global'),
                   # signed trasitivity (Kunegis, 2009)
                   s.trans = sum(net*(net%*%t(net)))/sum(abs(net)%*%t(abs(net)))
                   )

              # output 
              return(res)
          })

# function extracting genes in pathway
setGeneric('getGenes', function(object, ...) standardGeneric('getGenes'))

setMethod('getGenes', signature='pathway',
          definition = function(object){
              return(rownames(object@adj))
          })

# plot method
if (!isGeneric("plot")) setGeneric('plot')

setMethod('plot', signature(x='pathway',y='missing'),
          function(x, y=NA,
                   highlight.genes = NULL, gene.names=c('legend','nodes',NA), main=NULL, 
                   asp=0.95, vertex.size=11, vertex.color='khaki1', vertex.label.cex=0.8,
                   edge.width=2, edge.color='olivedrab4'){
                         
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
            main=main)
# add legend
if(gene.names == 'legend'){
  legend(,x=-1.55, y=1.5, bty='n',lty=NULL, legend=paste(1:vcount(g),V(g)$name),cex=0.7)
}

invisible(NULL)
})

# function randomly selecting effect genes in pathway
setGeneric('sample_genes', function(object, ...) standardGeneric('sample_genes'))

setMethod('sample_genes', signature='pathway',
          definition = function(object, no=3){
            g <- pathway2igraph(object)

# sample a gene with high betweennes centrality
sel1 <- sample(V(g), size=1, prob=betweenness(g))
# sample two of its neighbors
sel2 <- sample(neighborhood(g,1,sel1)[[1]], size=no-1)
# combine samples
samp <- c(as.numeric(sel1),sel2)

            return(samp)
          })


#------- functions for updating pathways ------------

# function to get information file for a pathway (= a list of all genes in the 
# pathway with start and end points
pathway_info <- function(id,  ...){
## id: character, pathway identifier, i.e. hsa00123 as in KEGG database

    #get genes from kegg
    input <- scan(url(paste("http://togows.dbcls.jp/entry/pathway/"
                  ,id,"/genes",sep="")), what="character")
    #get gene names only
    g <- lapply(input, function(x) if(substr(x,nchar(x),nchar(x))==";")
                {substr(x,1,nchar(x)-1)})
    g[sapply(g, is.null)] <- NULL
    
    #extract gene infos about pathway's genes from Ensembl
    ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    info <- getBM(attributes=c("start_position","end_position",
            "chromosome_name","hgnc_symbol"), filters=c("hgnc_symbol"),
            values=g, mart=ensembl)
    info$chromosome_name <- as.numeric(as.character(info$chromosome_name))
    info <- na.omit(info)
    pathway_info <- cbind(rep(paste(id,sep=""),length(info[,1])), info)
    colnames(pathway_info) <- c("pathway","gene_start","gene_end","chr",
                              "gene")
    return(pathway_info)
}



#2 helper functions used in 'get_network_matrix'
set_one <- function(M){
  if(length(M[M>1])>0){ 
     print(paste("Edges value > 1 set to 1!",sep=""))
     M[M>1] <- 1 }
  if(length(M[M < (-1)])>0){
     print(paste(M, ": Edges value < -1 set to -1!",sep=""))
     M[M < (-1)] <- -1 }
  return(M)
}

#translate numbers to gene names:
set_names <- function(N,nodes,liste){
  name <- substr(nodes,5,nchar(nodes)) 
  for(i in 1:length(name)){ 
    name[i] <- liste[liste[,1]==name[i],2] 
  }
  colnames(N) <- rownames(N) <- name
  return(N)
}
   
# function to calculate the network matrix          
get_network_matrix <- function(id, directed,  ...){
## hsa: pathway identifier as in KEGG database, character
## directed: Should networkmatrix be directed? TRUE/FALSE

    #download KGML from KEGG:
    retrieveKGML(substr(id,4,nchar(id)), organism="hsa",
                 destfile=paste(id,".xml",sep=""), method="internal")
    filename <- paste(id,".xml",sep="")
  
    #get genes names and their numbers in kegg
    info <- scan(url(paste("http://togows.dbcls.jp/entry/pathway/",
                           id,"/genes",sep="")), what="character")   
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

    pathgraph <- parseKGML2Graph(filename, expandGenes=TRUE)
    nodes <- nodes(pathgraph)  #Vektor mit Nummern der Gene (format: "hsa:226")
    
    #get edges via parsing to edgelist:
    edgelist <- parseKGML2DataFrame(filename)
    edgelist <- edgelist[!is.na(edgelist[,3]),] #delete NA-types            

  # --- Skip pathway, if no edges or wrong number of genes ---
    if(length(edgelist)==0){
       print(paste("KGML-file for ",id," has no edges!",sep=""))
       next}   
    if(length(nodes)!=length(liste[,2])){
       paste(paste("Wrong number of genes in ",id,"!",sep=""))
       next}               
    #set up empty matrix:
    N <- matrix(rep(0, (length(nodes)^2)),nrow=length(nodes))
    colnames(N) <- nodes
    rownames(N) <- nodes   
    verb.a <- edgelist[edgelist[,3]=="activation",] 
    verb.i <- edgelist[edgelist[,3]=="inhibition",]
    verb.s <- edgelist[edgelist[,3]!="inhibition"&edgelist[,3]!="activation",] 
  
  # --- if activations/inhibitions specified: Signed matrix ---
    if( (length(verb.a[,1])+length(verb.i[,1]))>0 ){
      if(length(verb.s[,1])==0){ print(paste(id,": Activation/Inhibition only: 
                                 Signed graph!",sep="")) }
      if(length(verb.s[,1])>0){ print(paste(id," has both: A/I edges and edges 
                                without type!",sep=""))}
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
       print(paste(hsa,": No edge-types, unsigned graph only!",sep=""))
       
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
      return(N)
    }else{
      return(M)
    }    
}
          
          
          
          
          
          

