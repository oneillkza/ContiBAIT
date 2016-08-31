locateMisorients.func <- function(compiledGrange, gapFile=NULL, stateNum=3, readCutOff=40, verbose=TRUE )
{

  compiledGrange <- split(compiledGrange, seqnames(compiledGrange), drop=TRUE)

  #only take gaps from chromosomes present in compiledGrange
  if(!(is.null(gapFile)))
  {
    gapFile <- gapFile[which(seqnames(gapFile) %in% names(compiledGrange))]
    gapFile <- split(gapFile, seqnames(gapFile), drop=TRUE)
  }

  fullMisoData <- GRanges()
  fullChimData <- GRanges()

  for(num in seq_len(length(names(compiledGrange))) )
  {
    ch <- names(compiledGrange)[num]
    if(verbose){message("Looking for issues on ", ch, " [",  which(names(compiledGrange) == ch), "/", length(compiledGrange), "]")}
    oneChr <- compiledGrange[[ch]]
    plusString <- as.character(strand(oneChr))
    
    plusString <- replace(plusString, plusString == "+", 1)
    plusString <- replace(plusString, plusString == "-", 0)
    plusString <- as.numeric(plusString)
    
    CNA.object <- suppressWarnings(CNA(plusString, as.character(seqnames(oneChr)), start(oneChr)))
    smoothed.CNA.object <- smooth.CNA(CNA.object, smooth.region=2)
    segmented <- segment(smoothed.CNA.object, verbose=0)
    
    segs <- segmented$output

  	rounder <- seq(0,1, length.out=stateNum)[2]	
  	pred.calls <- round(segs$seg.mean/rounder)*rounder
    rownames(segs) <- seq_len(nrow(segs))
    segs <- data.frame(chr=segs[,2], start=segs[,3], end=segs[,4], count=segs[,5], calls=pred.calls)
    
    f <- 2
    while (f < nrow(segs) )
    {
      #If fragment has too few reads to make accurate call in non-merged data, combine with upstream fragment
    	if(segs$count[f] <= readCutOff )
    	{
        #if region above is 1 and region below is 0 or visa versa: ie 0, 0.5, 1 or 1, 0.5, 0 or if 
    		if( (abs(segs$calls[f-1]-segs$calls[f+1]) == 1) | ((abs(segs$calls[f-1]-segs$calls[f+1]) != 0) & (segs$calls[f-1] == 0.5)) )
    		{
          #Make end of previous region to end of small region (ie combine these fragments)
          segs$end[f-1] <- segs$end[f]
          #And add to the count
          segs$count[f-1] <- sum(segs$count[(f-1):f])
          #then remove small fragment that has been merged
          segs <- segs[-f,]
    
     		}
        #If small region is flanked by the same number, eg a WC region between two WW regions, combine all three. ie 1,0,1 or 0,1,0 or 0.5,1,0.5 or 0.5,0,0.5 or 1,0.5,1 or 0,0.5,0
        else if(abs(segs$calls[f-1]-segs$calls[f+1]) == 0 )
        {
          #make start of upstream fragment to start of downstream (ie range now covers all three framents)
          segs$start[f+1] <- segs$start[f-1]
          #combine all three counts
          segs$count[f+1] <- sum(segs$count[(f-1):(f+1)])
          # and remove the two joined fragments
          segs <- segs[-c((f-1),f),]
     		}else{ 
          #if segs$calls[f+1] == 0.5, or if stateNum != 3 so these calculations don't work, combine with the downstream fragment
          segs$start[f+1] <- segs$start[f]
          #And add to the count
          segs$count[f+1] <- sum(segs$count[f:(f+1)])
          #then remove small fragment that has been merged
          segs <- segs[-f,]
        }

    	}else{
    		f <- f+1 
    	}
    }
    
    #If gapfile is seleced, use to extend regions and attempt span gaps between misorientations
    if(!(is.null(gapFile)))
    {
      chrGap <- gapFile[[ch]]
      if(!(is.null(chrGap)))
        {
        for(gap in seq(2:nrow(segs)))
        { 
          betweenGap <- chrGap[which(start(chrGap) >= segs$end[gap-1] & end(chrGap) <= segs$start[gap]),]
        	if(length(betweenGap) > 0)
        	{
        		if(segs$calls[(gap-1)] == 1)
        		{
        			#if first element is a positive call and second element is negative or WC, extend end position to end of last gap, and start position of next segment to end of last gap
        			segs$end[(gap-1)] <- end(betweenGap)[length(betweenGap)]
        			segs$start[gap] <- end(betweenGap)[length(betweenGap)]+1
        
        		}else if ( segs$calls[gap] == 1)
        		{
        			#if second element is a positive call, change start position to start of first gap
        			segs$start[gap] <- start(betweenGap)[1]
        			segs$end[(gap-1)] <- start(betweenGap)[1]-1
        		}
        	}
        }
      }
    }
    

    #remove 0.5 calls as these are probably not present on this chromosome, but are chimeric
    chimeraCalls <- segs[which(segs$calls == 0.5),]
    if(nrow(chimeraCalls) > 0)
    {
      chimeraCalls$name <- paste('chrUn_', chimeraCalls$chr, ':', chimeraCalls$start, '-', chimeraCalls$end , sep='')
      chimRanges <- GRanges(seqnames=chimeraCalls$chr, IRanges(start=chimeraCalls$start, end=chimeraCalls$end), name=chimeraCalls$name) #, score=chimeraCalls$count)
    }else{ 
      chimRanges <- GRanges()
    }

     chrCalls <- segs[which(segs$calls == 0),]
     if(nrow(chrCalls) > 0)
     {
        chrCalls$name <- paste(chrCalls$chr, ':', chrCalls$start, '-', chrCalls$end, sep="")

       	gCalls <- GRanges(seqnames=chrCalls$chr, IRanges(start=chrCalls$start, end=chrCalls$end), name=chrCalls$name) #, score=chrCalls$count)
        strand(gCalls) <- '-'

        fullMisoData <- suppressWarnings(append(fullMisoData, gCalls))
        fullChimData <- suppressWarnings(append(fullChimData, chimRanges))
      }
  }

  return(list(fullMisoData, fullChimData))
}  

####################################################################################################
#' locateMisorients -- function to identify libraries that hare similar WC patterns on chromosomes 
#' 
#' @param compiledGrange A GRanges object consisting of read locations. Can be an individual file or the product of thoroughBed
#' @param gapFile A GRanges object consisting of start and end locations of assembly gaps (defaul it NULL)
#' @param stateNum The number of expected strand states. Default is 3 (WW, WC and CC). Function may exhibit unusual behaviour is changed
#' @param readCutOff The minimal number of reads required to make an accurate strand state call. Default is 40.
#' @param verbose prints messages to the terminal (default is TRUE)
#' 
#' @return a list of ChrTable objects. The first is a ChrTable of misorientations detected. The second is a ChrTable of chimera detected.
#' Output can be used in downstream functions (strandSeqFreqTable)
#' @aliases locateMisorients locateMisorients,locateMisorients-GRanges-method
#' @rdname locateMisorients
#' @import Rsamtools
#' @import IRanges
#' @import GenomicRanges
#' @import DNAcopy
#' @importFrom S4Vectors DataFrame
#' @export
#' @include AllClasses.R
####################################################################################################

setMethod('locateMisorients',
          signature=signature(compiledGrange='GRanges'),
          definition = locateMisorients.func
          )
