plotContigOrder.func <- function(contigOrder, lg='all', verbose=TRUE)
{

	masterGroups <- sapply(1:nrow(contigOrder), function(x) strsplit(as.character(contigOrder[,1]), "\\.")[[x]][1])

	if(lg == 'all'){lg <- seq(1:length(unique(masterGroups)))}
	for(link in lg)
	{
		if(verbose){message(' -> Processing ', link)}
		contigOrderGrp <- contigOrder[grep(paste(unique(masterGroups)[link],"\\.", sep=""), contigOrder[,1]),]
#if(KeepZeros){rownames(contigOrderGrp) <- contigOrderGrp[,1]}

		if(nrow(as.matrix(contigOrderGrp)) > 2)
		{
			contigChr <- sub(':.*', '', contigOrderGrp[,2])
			primaryContigChr <- names(sort(table(contigChr), decreasing=TRUE))[1]

			contigLengths <- sub('.*:', '', contigOrderGrp[,2])
			contigStarts <- sub('-.*', '', contigLengths)

		if( length(unique(names(contigStarts))) != length(contigStarts))
		{
			#If more than one contig in the same sub-LG, take the mean start position.
			mergeFrame <- data.frame(lg=paste(names(contigStarts), contigChr, sep='LINK'), chr=contigChr, start=as.numeric(contigStarts)/10^6)
			mergeFrameAg <- aggregate(start~lg, mergeFrame, mean)
			contigOrderFrame <- mergeFrameAg[mergeFrame$lg,]
			contigOrderFrame <- data.frame(lg=sub('LINK.*', '', contigOrderFrame$lg), chr=sub('.*LINK', '', contigOrderFrame$lg), start=contigOrderFrame$start)
			contigOrderFrame$bin <- c(1:nrow(contigOrderFrame))
			contigOrderFrame$knownOrder <- as.numeric(rownames(contigOrderFrame[order(contigOrderFrame$start),]))
		}else{

			orderedLocation <- unlist(sapply(1:length(unique(contigOrderGrp[,1])), function(x) rep(x, length(contigOrderGrp[,1][which(contigOrderGrp[,1] == unique(contigOrderGrp[,1])[x])]))))
			contigOrderFrame <- data.frame(lg=names(contigChr), chr=contigChr, start=as.numeric(contigStarts)/10^6, bin=orderedLocation)
			contigOrderFrame$knownOrder <- as.numeric(rownames(contigOrderFrame[order(contigOrderFrame$start),]))
		}


			spearmanCor <- cor(contigOrderFrame$bin[which(contigOrderFrame$chr == primaryContigChr)], 
							   contigOrderFrame$knownOrder[which(contigOrderFrame$chr == primaryContigChr)], 
							   use="everything", 
							   method="spearman")

			if(spearmanCor < 0)
			{
				contigOrderFrame[,3:4] <- contigOrderFrame[nrow(contigOrderFrame):1, 3:4]
				spearmanCor <- spearmanCor*-1
			}		
	
#			rsquareOrient <- summary(lm(bin ~ knownOrder, data=contigOrderFrame[which(contigOrderFrame$chr == primaryContigChr),]))$r.squared
#			if(is.nan(rsquareOrient)){rsquareOrient <- 0}
#			rsquare <- summary(lm(bin ~ -1+knownOrder, data=contigOrderFrame[which(contigOrderFrame$chr == primaryContigChr),]))$r.squared
#			if(is.nan(rsquare)){rsquare <- 0}

#			if(rsquare < rsquareOrient)
#			{
#				contigOrderFrame[,3:4] <- contigOrderFrame[nrow(contigOrderFrame):1, 3:4]
#				rsquare <- rsquareOrient
#			}

			#rsquare <- round(rsquare, digits=2)
			spearmanCor <- round(spearmanCor, digits=2)
			print(ggplot(contigOrderFrame, aes_string("bin", "start") )+
			geom_point(aes_string(x="bin", y="start" , colour="chr"), size=2)+
			labs(x="contiBAIT predicted location of contigs", y="Assembly ordered location of contigs (Mb)")+
			geom_smooth(method="lm")+
			ggtitle(paste(primaryContigChr, 
						 " plot of ", 
						 length(contigChr), 
						 " fragments (", 
						 length(unique(contigOrderFrame$bin)), 
						 " sub-linkage groups)\nSpearman correlation = ", 
						 spearmanCor, 
						 sep="")))
		}
	}
}
####################################################################################################
#' Plot ordering of contigs within a single linkage group.
#' @param contigOrder matrix from orderAllContigs with the subdivided linkage groups and the names of the contigs to plot
#' @param lg Integer specifying the linkage group by which to plot. Default is all
#' @param verbose prints messages to the terminal (default is TRUE)
#' @aliases plotContigOrder plotContigOrder,ContigOrdering-method
#' @rdname plotContigOrder
#' @import ggplot2
#' @example inst/examples/plotContigOrder.R
#' @return A ggplot object (which will be plotted automatically if not assigned).
#' @export
####################################################################################################

setMethod('plotContigOrder',
      signature = signature(contigOrder='ContigOrdering'),
      definition = plotContigOrder.func
      )
