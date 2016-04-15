fixLinkageGroups.func <- function(contigOrdering, orderFrame, linkageGroupList, whichLG=NULL, relatedCutOff=0.6, verbose=TRUE)
{
	if(is.null(whichLG)){whichLG=seq_len(length(orderFrame))}

	#order the linkage groups first...
	if(verbose){message(" -> Ordering linkage groups")}
	for(lg in whichLG)
  	{
  		lgName <- paste("LG", lg, sep="")
  		orderMyLGs <- contigOrdering[which(rownames(contigOrdering) == lgName),2]
  		names(orderMyLGs) <- NULL
  		linkageGroupList[[lg]] <- orderMyLGs
  	}

  	newLGList <- list()
	rownames(contigOrdering) <- contigOrdering[,1]
    exdiag <- function(mat, off){mat[row(mat)+off == col(mat)]}

	for(lg in whichLG)
  	{
  		orderGroup <- orderFrame[[lg]]
  		chrNames <- rownames(orderGroup)
		lgOrdering <- contigOrdering[which(contigOrdering[,1] %in% rownames(orderGroup)),]
	    orderGroup <- data.frame(orderGroup)
	    orderGroup <-  data.frame(lapply(orderGroup, function(x) factor(x, levels=c(1,2,3))))
	    rownames(orderGroup) <- chrNames


		similarLinkageStrands <- as.matrix(1-daisy(orderGroup))
	    diag(similarLinkageStrands) <- 1
	    #sanity check that group seems fine...
	    diagOffset <- exdiag(similarLinkageStrands,1)
	 	breakpoints <- c(0, which(diagOffset < relatedCutOff), nrow(similarLinkageStrands))

	 	if(length(breakpoints) > 2)
	 	{
			if(verbose){message(" -> LG", lg, ": ", length(breakpoints)-1, " mis-clustered regions found")}

			for(num in seq_len(length(breakpoints)-1))
			{
				#for cases where there is a group of only one
				if((breakpoints[num]+1)-breakpoints[num+1] == 0)
				{
					splitStrands <- rownames(similarLinkageStrands)[breakpoints[num]+1]
					splitStrands <- lgOrdering[splitStrands,2]					
				}else{
					splitStrands <- similarLinkageStrands[(breakpoints[num]+1):breakpoints[num+1], (breakpoints[num]+1):breakpoints[num+1]]
					splitStrands <- lgOrdering[which(lgOrdering[,1] %in% rownames(splitStrands)),2]
					names(splitStrands) <- NULL
				}
				newLGList <- c(newLGList, list(splitStrands))
		    }
		}else{
			if(verbose){message(" -> LG", lg, ": Cluster verified")}
			newLGList <- c(newLGList, list(linkageGroupList[[lg]]))
		}
	}


	newLGList <- newLGList[order(sapply(newLGList, length), decreasing=TRUE)]
	#order linkage groups by biggest first
	newLGList <- LinkageGroupList(
  						 newLGList, 
  						  names= sapply(1:length(newLGList), 
  						  			  function(x){
  						  			  	paste('LG', x, ' (', length(newLGList[[x]]), ')', sep='') 
  						  			  	}))
	return(newLGList)
}

####################################################################################################
#' fixLinkageGroups -- searches for discrepancies within ordered contigs to highlight erroneously merged fragments.
#'
#' @param contigOrdering a data.frame of ordered contigs with linkage group names of class ContigOrdering
#' @param orderFrame a list of StrandStateMatrix elements of class StrandStateList
#' @param linkageGroupList List of vectors containing names of contigs belonging to each LG of type LinkageGroupList.
#' @param whichLG  vector of integers specifying the element(s) of linkageGroupList to be ordered (i.e. which specific linkage groups to try to order). Default is all LGs.
#' @param relatedCutOff The minimal dissimilarity between adjacent contigs to subset a linkage group into multiple smaller groups. Default is 0.6 
#' @param verbose Outputs information to the terminal. Default is TRUE.
#' @aliases fixLinkageGroups fixLinkageGroups-ContigOrdering-StrandStateList-method
#' @rdname fixLinkageGroups
#' @return a LinkageGroupList with erroneously clustered contigs seperated into their own groups
#' 
#' @importFrom cluster daisy
#' @export
#' @include AllClasses.R
####################################################################################################

 setMethod('fixLinkageGroups',
 		  signature = signature(contigOrdering='ContigOrdering',
 		  						orderFrame='StrandStateList',
 		  						linkageGroupList='LinkageGroupList'),
 		  definition = fixLinkageGroups.func
 		  )
