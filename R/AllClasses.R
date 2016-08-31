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
		 representation('list', names='character'),
		 validity=function(object){is.list(object) && (length(names(object)) >=1)}
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
#' A class for storing library group calls for contigs
#' 
#' \describe{
#'  This class is a list of lists, each primary list element is a chromosome/contig, and contains 2 sub-list elements: a list of 'Mostly Watson" and "Mostly Crick" library names. 
#' }
#
#' @export
#' @rdname LibraryGroupList

setClass("LibraryGroupList", 
		 representation('list', names='character'),
		 # Ensure all list elements in LibraryGroupList are of class LinkageGroupList
		 validity=function(object){unique(sapply(seq_along(object), function(x) class(object[[x]]) == "LinkageGroupList"))}
		 )

#' Constructor forLibraryGroupList
#' @aliases LibraryGroupList
#' @rdname LibraryGroupList
#' @param libraryGroups a list of lists, with each primary list element representing a chromosome with two internal list elements; a character vector of mostly watson library names, and a character vector of mostly Crick library names
#' @param names a vector of names of linkage groups
#' @return a \code{LibraryGroupList}
#' @export
#' @examples
#' lg1 <- LinkageGroupList(list(a=c('library1', 'library2'), b=c('library3')), names=c('chr1_Mostly_Crick', 'chr1_Mostly_Watson'))
#' lg2 <- LinkageGroupList(list(a=c('library1'), b=c('library6', 'library4')), names=c('chr2_Mostly_Crick', 'chr2_Mostly_Watson'))
#' libList <- LibraryGroupList(list(lg1, lg2))

LibraryGroupList <- function(libraryGroups = list(), names=character())
{
	names <- sapply(seq_along(libraryGroups), function(x) strsplit(names(libraryGroups[[x]][1]), "_")[[1]][1])	
	new('LibraryGroupList', libraryGroups, names=names)	
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
		 contains='GRanges')

#' Constructor for ChrTable
#' @aliases ChrTable
#' @rdname ChrTable
#' @param chrRanges  a GRanges object with a meta column 
#' called name, which represents the fragment name
#' @return a \code{ChrTable}
#' @export

ChrTable <- function(chrRanges=GRanges())
{
	if(is.null(chrRanges$name)){chrRanges$name <- paste(seqnames(chrRanges), ':', start(chrRanges), '-', end(chrRanges), sep='')}
	new('ChrTable', chrRanges)	
}

# =========================================================================
#' A class for storing StrandStateList lists for contigs
#' 
#' \describe{
#'  This class is a list of StrandStateMatrices, each a subset of the StrandStateMatrix split by linkagegroups
#' }
#'
#' @export
#' @rdname StrandStateList

setClass("StrandStateList", 
		 representation('list', names='character'),
		 # Ensure all list elements in LibraryGroupList are of class StrandStateMatrix
		 validity=function(object){unique(sapply(seq_along(object), function(x) class(object[[x]]) == "StrandStateMatrix"))}
		 )

#' Constructor StrandStateList
#' @aliases StrandStateList
#' @rdname StrandStateList
#' @param strandGroupList a list of StrandStateMatrix elements, with each primary element representing a StrandStateMatrix containing ordered contigs from a LinkageGroupList element
#' @param names a vector of names of StrandStateMatrix elements
#' @return a \code{StrandStateList}
#' @export

StrandStateList <- function(strandGroupList= list(), names=character())
{
	new('StrandStateList', strandGroupList, names=names)
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

#' Example of a LibraryGroupList, containing library names
#' @name exampleLibList
#' @docType data
#' @keywords data
NULL

#' Example of a ContigOrdering table, containing a list with a matrix of ordered groups element
#' and a StrandStateList element
#' @name exampleContigOrder
#' @docType data
#' @keywords data
NULL