## =========================================================================
## Generic for preprocessStrandTable.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export preprocessStrandTable
setGeneric("preprocessStrandTable", 
		   function(strandTable, 
		   			strandTableThreshold=NULL, 
		   			filterThreshold=NULL, 
		   			orderMethod=NULL, 
		   			lowQualThreshold=NULL, 
		   			verbose=NULL, 
		   			minLib=NULL, 
		   			ignoreInternalQual=NULL) standardGeneric("preprocessStrandTable"))


## =========================================================================
## Generic for clusterContigs.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export clusterContigs
setGeneric("clusterContigs", 
		   function(object, 
		   		 similarityCutoff=NULL,
		   		 recluster=NULL, 
		   		 minimumLibraryOverlap=NULL,
		   		 randomise=NULL,
		   		 randomSeed=NULL,
		   		 randomWeight=NULL,
		   		 snowCluster=NULL,
		   		 clusterBy=NULL,
		   		 verbose=NULL) standardGeneric("clusterContigs"))


## =========================================================================
## Generic for reorientLinkageGroups.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export reorientLinkageGroups
setGeneric("reorientLinkageGroups", 
		   function(object, 
		   			allStrands,
		   			previousOrient=NULL,
		   		 	verbose=NULL) standardGeneric("reorientLinkageGroups"))

## =========================================================================
## Generic for mergeLinkageGroups.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export mergeLinkageGroups
setGeneric("mergeLinkageGroups", 
           function(object, 
                    allStrands,
                    clusNum=NULL, 
                    cluster=NULL,
                    similarityCutoff=NULL) standardGeneric("mergeLinkageGroups"))

## =========================================================================
## Generic for orderAllLinkageGroups.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export orderAllLinkageGroups
setGeneric("orderAllLinkageGroups", 
           function(linkageGroupList,
            		strandStateMatrix, 
            		strandFreqMatrix, 
            		strandReadCount, 
            		whichLG=NULL, 
            		saveOrdered=NULL, 
            		orderCall=NULL, 
            		randomAttempts=NULL, 
            		verbose=NULL) standardGeneric("orderAllLinkageGroups"))

## =========================================================================
## Generic for plotLGDistances.R
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export plotLGDistances
setGeneric("plotLGDistances", 
           function(object, 
                    allStrands,
                    lg=NULL,
                    labels=NULL) standardGeneric("plotLGDistances"))

## =========================================================================
## Generic for barplotLinkageGroupCalls
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export barplotLinkageGroupCalls
setGeneric("barplotLinkageGroupCalls", 
           function(object, 
                    assemblyBED, 
                    by=NULL, 
                    returnTable=NULL) standardGeneric("barplotLinkageGroupCalls"))

## =========================================================================
## Generic for plotWCdistribution
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export plotWCdistribution
setGeneric("plotWCdistribution", 
           function(object, 
                    allStrands,
                    filterThreshold=NULL) standardGeneric("plotWCdistribution"))

## =========================================================================
## Generic for makeBoxPlot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export makeBoxPlot
setGeneric("makeBoxPlot", 
           function(chrTable, 
                    linkage.contigs) standardGeneric("makeBoxPlot"))

## =========================================================================
## Generic for ideogramPlot
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export ideogramPlot
setGeneric("ideogramPlot", 
           function(WatsonFreqList, 
                    CrickFreqList,
                    chrTable,
                    plotBy=NULL,
                    showPage=NULL,
                    orderFrame=NULL,
   		            orientationData=NULL,
                    verbose=NULL) standardGeneric("ideogramPlot"))

## =========================================================================
## Generic for writeBed
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' @export writeBed
setGeneric("writeBed",
		   function(chrTable, 
					orientationData, 
					contigOrder,
					libWeight=NULL,
					file=NULL) standardGeneric("writeBed"))

## =========================================================================
## show Methods
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

## show StrandStateMatrix
#' @name show,StrandStateMatrix-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a StrandStateMatrix
#' @return nothing
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
#' @return nothing
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
#' @return nothing
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
#' @return nothing
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
#' @return nothing
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
#' @return nothing
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

## show OrientationFrame
#' @name show,OrientationFrame-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a OrientationFrame
#' @return nothing
#' @description Shows a OrientationFrame
setMethod("show",
		  signature=signature(object="OrientationFrame"),
		  definition=function(object)
		  {
		  	elements <- nrow(object)
		  	misorientations <- nrow(object[which(object[,2] == '-'),])
		  	cat('A data.frame of ', elements, ' contigs with ',misorientations,' identified misorientations.\n')
		  }
)


## show ChrTable
#' @name show,ChrTable-method
#' @export
#' @docType methods
#' @title show-methods
#' @param object a ChrTable
#' @return nothing
#' @description Shows a ChrTable
setMethod("show",
		  signature=signature(object="ChrTable"),
		  definition=function(object)
		  {
		  	if(ncol(object) == 2)
		  	{
			  	cat('A data.frame of', length(object[,1]), 'fragments from a', sum(object[,2])/1000000, 'Mb genome.\n')
		  	}else{
		  		cat('A data.frame of', length(object[,1]), 'fragments from a', (sum(as.numeric(object[,3]))-sum(as.numeric(object[,2])))/1000000, 'Mb genome.\n')
		  	}
		  	show(head(object))
			cat('...\n')
			show(tail(object))
		  }
)