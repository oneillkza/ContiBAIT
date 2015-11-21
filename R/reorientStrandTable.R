

reorientStrandTable.func <- function(object, linkageGroups, orientation)
{
	reorientedStrands <- object
	
	toReorient <- unlist(linkageGroups[which(orientation== '-')])
	toReorientStrands <- suppressWarnings(reorientedStrands[toReorient,])
	
	toReorientStrands <- suppressWarnings(data.frame(
		apply(toReorientStrands, c(1,2), 
			  function(entry)
			  {
			  	if(!is.na(entry)&&entry=='3') return('1')
			  	if(!is.na(entry)&&entry=='1') return('3')
			  	entry
			  })
	))
	
	reorientedStrands <- data.frame(lapply(reorientedStrands, function(x){factor(x, levels=c(1,2,3))}))  
	rownames(reorientedStrands) <- rownames(object)
	
 	reorientedStrands[toReorient,] <- toReorientStrands
	suppressWarnings(return(new('StrandStateMatrix', reorientedStrands)))
}

####################################################################################################
#' Reorient a strand table, given the orientation vector computed from reorientLinkageGroups
#' @param object Strand table to reorient
#' @param linkageGroups List of vectors containing names of contigs belonging to each LG
#' @param orientation vector of linkage group orientations, as '+' or '-'
#' @return the strand table, reoriented according to the linkage group orientations
#' @aliases reorientStrandTable reorientStrandTable,StrandStateMatrix,StrandStateMatrix-method
#' @example inst/examples/clusterContigs.R
#' @export
####################################################################################################

setMethod('reorientStrandTable',
		  signature = signature(object='StrandStateMatrix'),
		  definition = reorientStrandTable.func
)