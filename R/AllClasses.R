# This file also contains roxygen entries for example data

# =========================================================================
#' A class for storing a matrix of read count-based strand states 
#' of a set of contigs over several libraries
#' 
#' \describe{
#'  The strand information stored in this object is the ratio of Watson to Crick reads
#'  mapping to each contig in each library (cell). This should fall within the range (-1,1).
#'	This class simply extends matrix, but with additional validity checking.
#' }
#
#' @export
#' @rdname StrandMatrix

setClass("StrandMatrix", 
		 contains='matrix', 
		 validity=function(object){is.double(object)})



StrandStateMatrixValidity <- function(object)
{
	areFactors <- all(sapply(object, function(x){is.factor(x)}))
  haveCorrectLevels <- all(sapply(object, function(this.col) {
     setequal(levels(this.col),as.factor(c('1','2','3')))
    }))
  areFactors && haveCorrectLevels
}

# =========================================================================
#' A class for storing a data frame of discrete strand states 
#' of a set of contigs over several libraries
#' 
#' \describe{
#'  The strand information stored in this object is a call of the strand state
#'  of each contig in each library.
#'  mapping to each contig in each library (cell). This should fall within the range (-1,1).
#'	This class simply extends matrix, but with additional validity checking.
#' }
#
#' @export
#' @rdname StrandStateMatrix


setClass("StrandStateMatrix", 
		 contains='data.frame', 
		 validity=StrandStateMatrixValidity)

# =========================================================================
#' A class for storing a data frame of directional reads from a single contig
#' of a single library
#' 
#' \describe{
#'  The strand information stored in this object represent every read meeting the threshold 
#'  critera from the specified contig for a given library.
#' }
#
#' @export
#' @rdname RawReadStrands


setClass("RawReadStrands", 
		 contains='data.frame', 
		 )


# =========================================================================
#' A class for storing linkage group calls for contigs
#' 
#' \describe{
#'  This class is simply a list of character strings containing the names of linkage groups. 
#' }
#
#' @export
#' @rdname LinkageGroupList

setClass("LinkageGroupList", 
		 contains='list')


#' An example StrandStateMatrix containing WW, CC and WC calls for contigs
#'
#' @name exampleWCMatrix
#' @docType data
# @author Kieran O'Neill \email{blahblah@@roxygen.org}
#' @keywords data
NULL

#' An example StrandStateMatrix containing only WW and CC calls for contigs
#'
#' @name exampleWWCCMatrix
#' @docType data
# @author Kieran O'Neill \email{blahblah@@roxygen.org}
#' @keywords data
NULL

# =========================================================================
#' A class for storing contig ordering of a linkage group
#'
#' \describe{
#'  This class is list of two character vectors that represent the calculated ordering of a linkage group. 
#'  The first element of this list is the names of each contig in the calculated order.
#'  The second element is the names of the merged contigs, i.e. those with equal strand states across all libraries,
#'  in the calculated order.
#' }
#
#' @export
#' @rdname ContigOrdering

setClass("ContigOrdering",
         contains='list')