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

# show method
setMethod('show', signature='pathway',
          definition = function(object){
	      # summarize pathway information
              cat('An object of class ', class(object), ' with id ',object@id,'\n\n',sep='')
              cat('Pathway with adjacency matrix of dimension ', dim(object@adj)[1], ': \n',sep='')
	      # print pathway adjacency matrix
	      N <- object@adj
	      N[N!=0] <- object@sign
              print(N)
              invisible(NULL)
          })

# summary method
setGeneric('summary', function(object, ...) standardGeneric('summary'))

setMethod('summary', signature='pathway',
          definition = function(object){
              cat('An object of class ', class(object), '\n\n',sep='')
              cat('?? nodes and ?? links; ?? activations and ?? inhibitions \n',sep='')
              cat('Average degree: \n')
              invisible(NULL)
          })

# analyze method - analyze pathway network properties
setGeneric('analyze', function(object, ...) standardGeneric('analyze'))

setMethod('analyze', signature='pathway',
          definition = function(object){
              cat('An object of class ', class(object), '\n\n',sep='')
              cat('?? nodes and ?? links; ?? activations and ?? inhibitions \n',sep='')
              cat('Average degree: \n')
              invisible(NULL)
          })

# plot method
setGeneric('plot', function(object, ...) standardGeneric('plot'))

setMethod('plot', signature='pathway',
          definition = function(object, ...){
             print('define igraph plot method!!')
              invisible(NULL)
          })


