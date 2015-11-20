## =========================================================================
## Generic for clusterContigs
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export clusterContigs
setGeneric("clusterContigs", 
		   function(object, similarityCutoff=NULL,
		   		 recluster=NULL, 
		   		 minimumLibraryOverlap=NULL,
		   		 randomise=NULL,
		   		 randomSeed=NULL,
		   		 randomWeight=NULL,
		   		 snowCluster=NULL,
		   		 verbose=NULL) standardGeneric("clusterContigs"))


## =========================================================================
## Generic for reorientStrandTable
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export reorientStrandTable
setGeneric("reorientStrandTable", 
		   function(object, 
		   		 linkageGroups=NULL, 
		   		 orientation=NULL) standardGeneric("reorientStrandTable"))

## =========================================================================
## Generic for mergeLinkageGroups.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export mergeLinkageGroups
setGeneric("mergeLinkageGroups", 
           function(object, 
                    allStrands, 
                    similarityCutoff=NULL) standardGeneric("mergeLinkageGroups"))

## =========================================================================
## Generic for plotLGDistances.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export plotLGDistances
setGeneric("plotLGDistances", 
           function(object, 
                    allStrands) standardGeneric("plotLGDistances"))

## =========================================================================
## show Methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## show StrandStateMatrix
#' @name show,StrandStateMatrix-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a StrandStateMatrix
setMethod("show",
		  signature=signature(object="StrandStateMatrix"),
		  definition=function(object)
		  {
        d <- dim(object)
  
		  	cat('A strand state matrix for ', d[[1]], ' contigs over ',d[[2]],' libraries.\n')
		  }
)

## show LinkageGroupList
#' @name show,LinkageGroupList-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a LinkageGrouplist
setMethod("show",
		  signature=signature(object="LinkageGroupList"),
		  definition=function(object)
		  {
		  	cat('A linkage group list containing ', length(object), ' linkage groups.\n\n')
		  	show(data.frame(NumberOfContigs=sapply(object, length)))
		  }
)

## show RawReadStrands
#' @name show,RawReadStrands-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a RawReadStrands
setMethod("show",
		  signature=signature(object="RawReadStrands"),
		  definition=function(object)
		  {
		  	cat('A data.frame containing ', nrow(object), ' reads.\n\n')
		  	#show(data.frame(NumberOfContigs=sapply(object, length)))
		  }
)