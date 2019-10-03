# Copyright (c) 2015, Mark Hills
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#    * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


####################################################################################################
#' mapGapFromOverlap -- function to co-localize strand state changes with assembly gaps 
#' 
#' @param sceFile GRanges object of strand state change locations in BED format 
#' @param gapFile GRanges object of assembly gaps in BED format (can be downloaded from UCSC table browser)
#' @param chrTable GRanges object of chromosome table (product of makeChrTable)
#' @param verbose prints messages to the terminal (default is TRUE)
#' @param overlapNum  Minimal number of strand state changes that overlap with a gap before assembly is cut at that location
#' 
#' @return a GRanges object of all contigs split by regions where the sceFile and gapFile GRanges objects overlap.
#' @import GenomicRanges
#' @import DNAcopy
#' @importFrom BiocGenerics "%in%"
#' @importFrom S4Vectors DataFrame
#' @export
#' @include AllClasses.R
####################################################################################################


mapGapFromOverlap <- function(sceFile,  gapFile, chrTable, verbose=TRUE, overlapNum=4)
{

	#remove names element from chrTable if present
	chrTable$name <- NULL

	hits <- countOverlaps(gapFile, sceFile)

	#append hits to gapFile to show number of SCE events that are coincident with gaps
	gapFile$overlapSCE <- hits

	#if not already, order by chromosome then start location of gap.
	gapFile <- sort(gapFile)

	# For every chromosome that has a gap and a corresponding SCE...
	for(chr in unique(seqnames(sceFile)[seqnames(sceFile) %in% seqnames(gapFile)])  )
	{
		gapPerChr <- gapFile[seqnames(gapFile) == as.character(chr)]
		if(max(gapPerChr$overlapSCE) >= overlapNum)
		{
			CNA.object <- CNA(gapPerChr$overlapSCE, 
							  as.character(seqnames(gapPerChr)), start(gapPerChr))
			smoothed.CNA.object <- smooth.CNA(CNA.object, smooth.region=2)
			segmented <- segment(smoothed.CNA.object, verbose=0)
			segs <- segmented$output
			if(nrow(segs) > 1)
			{
				segs <- segs[which(diff(sign(diff(segs$seg.mean)))==-2)+1,c(2:4, 6)]
				segs <- segs[which(segs$seg.mean >= overlapNum),]
				segs <- GRanges(segs$chrom, IRanges(start=segs$loc.start, end=segs$loc.end))
			}else{
				segs <- gapPerChr[gapPerChr$overlapSCE == max(gapPerChr$overlapSCE)]
			}

			segs$overlapSCE <- NULL
			chrTable <- append(GRanges(chrTable), segs)
		}
	}

	#Now split the assembly where overlaps occur
	chrTable <- disjoin(chrTable)

	chrTable$name <- paste(seqnames(chrTable), ":", start(chrTable), "-", end(chrTable), sep="")
	chrTable <- ChrTable(chrTable)
	return(chrTable)
}
