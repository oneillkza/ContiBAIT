####################################################################################################
#' Plot ordering of contigs within a single linkage group.
#' @param contigOrder data.frame from orderAllContigs with the subdivided linkage groups and the names of the contigs to plot
#' @param lg Integer specifying the linkage group by which to plot
#' @import ggplot2
#' @example inst/examples/plotContigOrder.R
#' @export
####################################################################################################

plotContigOrder <- function(contigOrder, lg)
{

	contigOrder <- contigOrder[grep(paste(unique(sapply(1:nrow(contigOrder), function(x) strsplit(as.character(contigOrder$LG), "\\.")[[x]][1]))[lg],"\\.", sep="") , contigOrder$LG),]

	contigChr <- sub(':.*', '', contigOrder[,2])
	primaryContigChr <- names(sort(table(contigChr), decreasing=T))[1]

	contigLengths <- sub('.*:', '', contigOrder[,2])
	contigStarts <- sub('-.*', '', contigLengths)

	orderedLocation <- unlist(sapply(1:length(unique(contigOrder[,1])), function(x) rep(x, length(contigOrder[,1][which(contigOrder[,1] == unique(contigOrder[,1])[x])]))))

	contigOrderFrame <- data.frame(chr=contigChr, start=as.numeric(contigStarts)/10^6, bin=orderedLocation)


	rsquareOrient <- summary(lm(bin ~ start, data=contigOrderFrame))$r.squared
	rsquare <- summary(lm(bin ~ -1+start, data=contigOrderFrame))$r.squared
	if(rsquare < rsquareOrient)
	{
		contigOrderFrame[,1:2] <- contigOrderFrame[nrow(contigOrderFrame):1, 1:2]
		rsquare <- summary(lm(bin ~ 1+start, data=contigOrderFrame))$r.squared
	}

	rsquare <- round(rsquare, digits=2)
	ggplot(contigOrderFrame, aes(bin, start) )+
	geom_point(aes(x=bin, y=start , colour=chr), size=2)+
	labs(x="contiBAIT predicted location of contigs", y="Assembly ordered location of contigs (Mb)")+
	geom_smooth(method="lm", se=FALSE, formula=y~x-1)+
	ggtitle(paste("Plot of ", length(contigChr), " fragments (", length(unique(contigOrderFrame$bin)), " sub-linkage groups)\nR-squared = ", rsquare,  sep=""))
}