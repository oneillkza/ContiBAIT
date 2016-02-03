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
#' @import methods
#' @import grDevices

setClass("StrandReadMatrix", 
		 contains='matrix', 
		 validity=function(object){is.integer(object)})

#' Constructor for StrandReadMatrix
#' @aliases StrandReadMatrix
#' @rdname StrandReadMatrix
#' @param counts an integer matrix of read counts
#' @return a \code{StrandReadMatrix}
#' @export
#' @examples
#' data("exampleWatsonFreq")
#' StrandReadMatrix(exampleWatsonFreq[,2, drop=FALSE])

StrandReadMatrix <- function(counts = matrix(integer()))
{
	new('StrandReadMatrix', counts)	
}


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


#' Constructor for StrandFreqMatrix
#' @aliases StrandFreqMatrix
#' @rdname StrandFreqMatrix
#' @param counts a double matrix of read count ratios
#' @return a \code{StrandFreqMatrix}
#' @export
#' @examples
#' data("exampleWatsonFreq")
#' data("exampleCrickFreq")
#' frequencyMatrix <- sapply(1:ncol(exampleCrickFreq), 
#' function(colNum){exampleCrickFreq[,colNum] / exampleWatsonFreq[,colNum]})
#' 
#' StrandFreqMatrix(frequencyMatrix)

StrandFreqMatrix <- function(counts = matrix(double()))
{
	new('StrandFreqMatrix', counts)	
}


StrandStateMatrixValidity <- function(object)
{
	obVec <- as.vector(object)
	dataCount <- length(grep("[123]", obVec))
	naCount <- length(which(is.na(obVec)))
	dataCount+naCount == length(obVec)
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
		 contains='matrix', 
		 validity=StrandStateMatrixValidity
		 )


#' Constructor for StrandStateMatrix
#' @aliases StrandStateMatrix
#' @rdname StrandStateMatrix
#' @param states an integer matrix of strand states by library
#' @return a \code{StrandStateMatrix}
#' @export
#' @examples
#' StrandStateMatrix(matrix(ncol=2, c(1,3,1,2)))

StrandStateMatrix <- function(states = matrix(integer()))
{
	new('StrandStateMatrix', states)	
}



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

#' Constructor forLinkageGroupList
#' @aliases LinkageGroupList
#' @rdname LinkageGroupList
#' @param linkageGroups a list of character vectors of names of contigs in each LG
#' @param names a vector of names of linkage groups
#' @return a \code{LinkageGroupList}
#' @export
#' @examples
#' lgList <- LinkageGroupList(list(lg1=c('contig1', 'contig2'), lg2=c('contig3')),
#' 								names=c('lg1', 'lg20'))

LinkageGroupList <- function(linkageGroups = list(), names=character())
{
	new('LinkageGroupList', linkageGroups, names=names)	
}



# =========================================================================
#' A class for storing contig ordering of a linkage group
#'
#' \describe{
#'  This class is a matrix of two character vectors that represent the calculated ordering of a linkage group. 
#'  The first element of this matrix is the Linkage Group sub-setted by contigs with equal strand states across all libraries
#'  in the calculated order. 
#'  The second element is the names of names of each contig in the calculated order.
#' }
#
#' @export
#' @rdname ContigOrdering

setClass("ContigOrdering",
		 contains='matrix')


#' Constructor for ContigOrdering
#' @aliases ContigOrdering
#' @rdname ContigOrdering
#' @param ordering a matrix of two character vectors that represent the calculated ordering of a linkage group. 
#'  The first element of this matrix is the Linkage Group sub-setted by contigs with equal strand states across all libraries
#'  in the calculated order. 
#'  The second element is the names of names of each contig in the calculated order.
#' @return a \code{ContigOrdering}
#' @export
#' @examples
#' thisOrdering <- ContigOrdering(matrix(ncol=2, c( "LG1.11", "chr2:1000820-2001640", 
#' 			"LG1.1", "chr2:3002461-4003281")))

ContigOrdering <- function(ordering=character())
{
	new('ContigOrdering', ordering)	
}


# =========================================================================
#' A class for storing contig orientations
#'
#' \describe{
#'  This class is a matrix of two character vectors that represent the orientation of contigs. 
#'  The first element of thismatrix is the contigs name
#'  The second element is the orinetation (as either + or -).
#' }
#
#' @export
#' @rdname OrientationFrame

setClass("OrientationFrame",
		 contains='matrix')

#' Constructor for OrientationFrame
#' @aliases OrientationFrame
#' @rdname OrientationFrame
#' @param orientation  a matrix of two character vectors that represent the orientation of contigs. 
#'  The first element of thismatrix is the contigs name
#'  The second element is the orinetation (as either + or -).
#' @return a \code{OrientationFrame}
#' @export
#' @examples
#' OrientationFrame(matrix(ncol=2, c("chr4:3002423-4003230", "+", 
#' 		"chr4:140113083-141113889", "+")))

OrientationFrame <- function(orientation=character())
{
	new('OrientationFrame', orientation)	
}

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

#' Constructor for ChrTable
#' @aliases ChrTable
#' @rdname ChrTable
#' @param chrRanges  a GRanges object with a meta column 
#' called name, which represents the fragment name
#' @return a \code{ChrTable}
#' @export

ChrTable <- function(chrRanges=GRanges())
{
	new('ChrTable', chrRanges)	
}


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