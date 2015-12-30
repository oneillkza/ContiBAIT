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
#' @param chrTable data..frame of chromosome table (product of makeChrTable)
#' @param verbose prints messages to the terminal (default is TRUE)
#' @param overlapNum  Minimal number of strand state changes that overlap with a gap before assembly is cut at that location
#' 
#' @return a data,frame in bed format of all contigs split by regions where the sceFile and gapFile GRanges objects overlap.
#' @import GenomicRanges
#' @import DNAcopy 
#' @importFrom S4Vectors DataFrame
#' @export
#' @include AllClasses.R
####################################################################################################


mapGapFromOverlap <- function(sceFile,  gapFile, chrTable, verbose=TRUE, overlapNum=4)
{
	if(ncol(chrTable) == 2 )
	{
		chrTable[,3] <- chrTable[,2]
		chrTable[,2] <- 0
	}
#	gapFileGRange <- GRanges(gapFile[,1], IRanges(gapFile[,2], gapFile[,3]))
#	sceFileGRange <- GRanges(sceFile[,1], IRanges(sceFile[,2], sceFile[,3]))
	hits <- countOverlaps(gapFile, sceFile)

	#append hits to gapFile to show number of SCE events that are coincident with gaps
	gapFile <- as.data.frame(gapFile)
	gapFile$overlapSCE <- hits

	#if not already, order by chromosome then start location of gap.
	gapFile <- gapFile[order(gapFile[,1], gapFile[,2]),]
	gapOverlap <- data.frame(chr=vector(), start=vector(), end=vector() )
	# For every chromosome that has a gap and a corresponding SCE...
	for(chr in unique(gapFile[which(gapFile[,1] %in% as.character(seqnames(sceFile))),1]) )
	{
		gapPerChr <- gapFile[which(gapFile[,1] == chr),]
		if(max(gapPerChr$overlapSCE) >= overlapNum)
		{
			CNA.object <- CNA(gapPerChr$overlapSCE, gapPerChr[,1], gapPerChr[,2])
			smoothed.CNA.object <- smooth.CNA(CNA.object, smooth.region=2)
			segmented <- segment(smoothed.CNA.object, verbose=0)
			segs <- segmented$output
			if(nrow(segs) > 1)
			{
				segs <- segs[which(diff(sign(diff(segs$seg.mean)))==-2)+1,c(2:4, 6)]
				segs <- segs[which(segs$seg.mean >= overlapNum),]
				completeSegs <- data.frame(chr=c(chrTable[chrTable[,1] == chr,1], as.character(segs$chrom)), start=c(chrTable[chrTable[,1] == chr,2], segs$loc.end), end=c(segs$loc.start, chrTable[chrTable[,1] == chr,3]))
			}else{
				segs <- gapPerChr[which(gapPerChr$overlapSCE == max(gapPerChr$overlapSCE)),]
				completeSegs <- data.frame(chr=c(chrTable[chrTable[,1] == chr,1], as.character(segs$seqnames)), start=c(chrTable[chrTable[,1] == chr,2], segs$start), end=c(segs$end, chrTable[chrTable[,1] == chr,3]))
			}
			gapOverlap <- rbind(gapOverlap, completeSegs)
		}
	}

	#THESE ARE THE GAPS!!!
	# now, add back the chromosomes which haven't been split 
	toInclude <- chrTable[which(!(chrTable$chr %in% gapOverlap[,1])),]
	colnames(toInclude) <- colnames(gapOverlap)
	splitContigData <- rbind(gapOverlap, toInclude)
	#splitContigData <- cbind(splitContigData, contigName=paste(splitContigData$chr, ":", splitContigData$start, "-", splitContigData$end, sep=""))
	rownames(splitContigData) <- paste(splitContigData$chr, ":", splitContigData$start, "-", splitContigData$end, sep="")
	return(new('ChrTable', splitContigData))
}