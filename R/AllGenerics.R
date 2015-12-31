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
## Generic for barplotLinkageGroupCalls
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export barplotLinkageGroupCalls
setGeneric("barplotLinkageGroupCalls", 
           function(object, 
                    assemblyBED, 
                    by=NULL, 
                    returnTable=NULL, 
                    saveFile=NULL,
                    ...) standardGeneric("barplotLinkageGroupCalls"))

## =========================================================================
## Generic for plotWCdistribution
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export plotWCdistribution
setGeneric("plotWCdistribution", 
           function(object, 
                    allStrands,
                    filterThreshold=NULL,
                    saveFile=NULL) standardGeneric("plotWCdistribution"))

## =========================================================================
## show Methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## show StrandStateMatrix
#' @name show,StrandStateMatrix-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a StrandStateMatrix
#' @description Shows a StrandStateMatrix
setMethod("show",
		  signature=signature(object="StrandStateMatrix"),
		  definition=function(object)
		  {
        d <- dim(object)
  
		  	cat('A strand state matrix for ', d[[1]], ' contigs over ',d[[2]],' libraries.\n')
		  }
)


## show StrandFreqMatrix
#' @name show,StrandFreqMatrix-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a StrandFreqMatrix
#' @description Shows a StrandFreqMatrix
setMethod("show",
		  signature=signature(object="StrandFreqMatrix"),
		  definition=function(object)
		  {
		  	d <- dim(object)
		  	
		  	cat('A matrix of strand frequencies for ', d[[1]], ' contigs over ',d[[2]],' libraries.\n')
		  }
)

## show StrandReadMatrix
#' @name show,StrandReadMatrix-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a StrandReadMatrix
#' @description Shows a StrandReadMatrix
setMethod("show",
		  signature=signature(object="StrandReadMatrix"),
		  definition=function(object)
		  {
		  	d <- dim(object)
		  	
		  	cat('A matrix of read counts for ', d[[1]], ' contigs over ',d[[2]],' libraries.\n')
		  }
)


## show LinkageGroupList
#' @name show,LinkageGroupList-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a LinkageGrouplist
#' @description Shows a LinkageGroupList
setMethod("show",
		  signature=signature(object="LinkageGroupList"),
		  definition=function(object)
		  {
		  	cat('A linkage group list containing ', length(object), ' linkage groups.\n\n')
		  	show(data.frame(NumberOfContigs=head(sapply(object, length)), row.names=NULL))
		  	if(length(object) > 5)
            show(data.frame("...           "=tail(sapply(object, length)), row.names=seq(length(object)-5, length(object) )))

		  }
)

## show RawReadStrands
#' @name show,RawReadStrands-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a RawReadStrands
#' @description Shows a RawReadStrands
setMethod("show",
		  signature=signature(object="RawReadStrands"),
		  definition=function(object)
		  {
		  	cat('A data.frame containing ', nrow(object), ' reads.\n\n')
		  	#show(data.frame(NumberOfContigs=sapply(object, length)))
		  }
)

## show ContigOrdering
#' @name show,ContigOrdering-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a ContigOrdering
#' @description Shows a ContigOrdering
setMethod("show",
		  signature=signature(object="ContigOrdering"),
		  definition=function(object)
		  {

		  	cat('A data.frame of', length(unique(sapply(1:nrow(object), function(x) strsplit(as.character(object$LG), "\\.")[[x]][1]))), 'LGs split into', length(unique(sapply(1:nrow(object), function(x) strsplit(as.character(object$LG), "\\.")[[x]][2]))), 'sub-groups from', nrow(object), 'ordered fragments.\n')
		  	if(length(unique(sapply(1:nrow(object), function(x) strsplit(as.character(object$LG), "\\.")[[x]][2]))) > 25)
		  	{
			  	show(head(table(object[,1])))
			  	cat('...')
			  	show(tail(table(object[,1])))
		  	}else{
		  		show(table(object[,1]))
		  	}
		  }
)


## show ChrTable
#' @name show,ChrTable-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a ChrTable
#' @description Shows a ChrTable
setMethod("show",
		  signature=signature(object="ChrTable"),
		  definition=function(object)
		  {
		  	if(ncol(object) == 2)
		  	{
			  	cat('A data.frame of', length(object[,1]), 'fragments from a', sum(object[,2])/1000000, 'Mb genome.\n')
		  	}else{
		  		cat('A data.frame of', length(object[,1]), 'fragments from a', (sum(object[,3])-sum(object[,2]))/1000000, 'Mb genome.\n')
		  	}
		  	show(head(object))
			cat('...\n')
			show(tail(object))
		  }
)