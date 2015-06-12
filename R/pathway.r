#################################################
#
# pathway object functions
#
#################################################

# !! define update function from make.BioPax.r !!
# return(pathway(...))

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
pathway2igraph <- function(object){

    # define adjacency matrix
    net <- object@adj
    net[net!=0] <- object@sign

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
setGeneric('sampleGenes', function(object, ...) standardGeneric('sampleGenes'))

setMethod('sampleGenes', signature='pathway',
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

