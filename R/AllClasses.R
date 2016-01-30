# This file also contains roxygen entries for example data


# =========================================================================
#' A class for storing read counts for a set of contigs over several libraries
#' 
#' \describe{
#'  The information stored in this class is simple read counts, so should be integers >=0.
#' }
#
#' @export
#' @rdname StrandReadMatrix

setClass("StrandReadMatrix", 
		 contains='matrix', 
		 validity=function(object){is.integer(object)})



# =========================================================================
#' A class for storing a matrix of frequencies of Watson to Crick reads
#' for a set of contigs over several libraries
#' 
#' \describe{
#'  The strand information stored in this object is the ratio of Watson to Crick reads
#'  mapping to each contig in each library (cell). This should fall within the range (-1,1).
#'	This class simply extends matrix, but with additional validity checking.
#' }
#
#' @export
#' @rdname StrandFreqMatrix

setClass("StrandFreqMatrix", 
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
		 representation('list', names='character')
		 )



# =========================================================================
#' A class for storing contig ordering of a linkage group
#'
#' \describe{
#'  This class is data.frame of two character vectors that represent the calculated ordering of a linkage group. 
#'  The first element of this data.frame is the Linkage Group sub-setted by contigs with equal strand states across all libraries
#'  in the calculated order. 
#'  The second element is the names of names of each contig in the calculated order.
#' }
#
#' @export
#' @rdname ContigOrdering

setClass("ContigOrdering",
		 contains='data.frame')

# =========================================================================
#' A class for storing contig orientations
#'
#' \describe{
#'  This class is data.frame of two character vectors that represent the orientation of contigs. 
#'  The first element of this data.frame is the contigs name
#'  The second element is the orinetation (as either + or -.
#' }
#
#' @export
#' @rdname OrientationFrame

setClass("OrientationFrame",
		 contains='data.frame')


# =========================================================================
#' A class for storing chromosome/fragment lengths
#'
#' \describe{
#'  This class is a GRanges object with a meta column called name, which represents the fragment name.
#' }
#
#' @export
#' @rdname ChrTable

setClass("ChrTable",
		 contains='GRanges',
		 validity=function(object){length(object$name) > 0})



# ========================================================================
# Data sets
#
# ========================================================================

#' Example of strand frequencies extracted from BAMS by strandSeqFreqTable
#'
#' @name exampleStrandFreq
#' @docType data
#' @keywords data
NULL

#' Example of read counts extracted from BAMS by strandSeqFreqTable
#'
#' @name exampleReadCounts
#' @docType data
#' @keywords data
NULL

#' An example StrandStateMatrix containing WW, CC and WC calls for contigs
#'
#' @name exampleWCMatrix
#' @docType data
#' @keywords data
NULL

#' An example Crick strand frequences extracted from BAMS by strandSeqFreqTable where BAITtables=TRUE
#'
#' @name exampleCrickFreq
#' @docType data
#' @keywords data
NULL

#' An example Watson strand frequences extracted from BAMS by strandSeqFreqTable where BAITtables=TRUE
#'
#' @name exampleWatsonFreq
#' @docType data
#' @keywords data
NULL

#' Example of a LinkageGroupList output from clusterContigs
#'
#' @name exampleLGList
#' @docType data
#' @keywords data
NULL

#' Example of a ChromosomeTable, containing contigs and their lengths
#' @name exampleChrTable
#' @docType data
#' @keywords data
NULL

#' Example of a divided chromosome, containing contigs and their lengths
#' @name exampleDividedChr
#' @docType data
#' @keywords data
NULL